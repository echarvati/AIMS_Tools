import os, shutil
from .gmx import GmxSimulation
from ...wrapper.gmx import *
from ...analyzer.acf import get_std_out
from ...panedr import edr_to_df

class Nvt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'nvt'
        self.logs = []

    def build(self, ppf=None, minimize=False):
        pass

    def prepare(self, prior_job_dir=None, gro='npt.gro', top='topol.top', T=298, jobname=None,
                dt=0.001, nst_eq=int(4E5), nst_run=int(5E5), random_seed=-1, nst_edr=5, nst_trr=50,
                tcoupl='v-rescale', diff_gk=False, mstools_dir=None, **kwargs) -> [str]:
        if prior_job_dir is None:
            raise Exception('prior_job_dir is needed for NVT simulation')

        # Copy gro and topology files from prior NPT simulation
        shutil.copy(os.path.join(prior_job_dir, gro), '.')
        shutil.copy(os.path.join(prior_job_dir, top), '.')
        for f in os.listdir(prior_job_dir):
            if f.endswith('.itp'):
                shutil.copy(os.path.join(prior_job_dir, f), '.')
        # Scale gro box for NVT simulation
        # TODO the equilibration of NPT simulation is not considered here
        box = self.gmx.get_box(os.path.join(prior_job_dir, 'npt.edr'))
        self.gmx.scale_box(gro, gro, box)

        nprocs = self.jobmanager.nprocs
        commands = []

        # NVT equilibrium with Langevin thermostat and Berendsen barostat
        self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-eq.mdp', T=T, gen_seed=random_seed,
                                           nsteps=nst_eq, nstxtcout=0, pcoupl='berendsen')
        cmd = self.gmx.grompp(mdp='grompp-eq.mdp', gro=gro, top=top, tpr_out='eq.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='eq', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NVT production with Langevin thermostat and Parrinello-Rahman barostat
        self.gmx.prepare_mdp_from_template('t_nvt.mdp', mdp_out='grompp-nvt.mdp', T=T, nsteps=nst_run, tcoupl=tcoupl,
                                           nstenergy=nst_edr, nstxout=nst_trr, nstvout=nst_trr, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-nvt.mdp', gro='eq.gro', top=top, tpr_out='nvt.tpr',
                              cpt='eq.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='nvt', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        if diff_gk:
            # diffusion constant, do not used, very slow
            commands.append(self.gmx.trjconv('nvt.tpr', 'nvt.trr', 'traj.gro', end=500, get_cmd=True))
            commands.append(os.path.join(mstools_dir, 'mstools', 'cpp', 'diff-gk') + ' traj.gro')
        # viscosity
        commands.append(self.gmx.energy('nvt.edr', properties=['Pres-XY', 'Pres-XZ', 'Pres-YZ'], out='pressure.xvg', get_cmd=True))
        weight = 0.00
        temperature = self.gmx.get_temperature_from_mdp('grompp-nvt.mdp')
        volume = self.gmx.get_volume_from_gro('npt.gro')
        commands.append(
            os.path.join(mstools_dir, 'mstools', 'cpp', 'vis-gk') + ' pressure.xvg' + ' %f' % (volume) + ' %f' % (
                temperature) + ' %.2f' % (weight))

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands


    # analyze diffusion constant
    def analyze_diff(self, charge_list, n_mol_list, diff_gk=False):
        # get temperature and volume
        volume_and_stderr = [self.gmx.get_volume_from_gro('nvt.gro'), 0.]
        [temperature_and_stderr] = self.gmx.get_properties_stderr('nvt.edr', ['Temperature'])

        # calculate diffusion constant using Einstein relation
        diff_e_dict = {'System': get_std_out(list(self.gmx.diffusion('nvt.xtc', 'nvt.tpr')))}
        for i in range(len(n_mol_list)):
            mol_name = 'MO%i' % (i)
            diff_e_dict.update({mol_name: get_std_out(list(self.gmx.diffusion('nvt.xtc', 'nvt.tpr', group=mol_name)))})

        info_dict = {'diffusion constant': diff_e_dict}

        # estimate electrical conductivity using Nernst-Einstein relation
        if charge_list != None and set(charge_list) != {0}:
            econ = 0.
            econ_stderr = 0.
            for i, charge in enumerate(charge_list):
                mol_name = 'MO%i' % (i)
                diff, stderr = diff_e_dict.get(mol_name)
                econ += diff * charge_list[i]**2 * n_mol_list[i]
                econ_stderr += stderr * charge_list[i]**2 * n_mol_list[i]
            econ *= 1.6 ** 2 / 1.38 * 10 ** 8 / temperature_and_stderr[0] / volume_and_stderr[0]
            econ_stderr *= 1.6 ** 2 / 1.38 * 10 ** 8 / temperature_and_stderr[0] / volume_and_stderr[0]
            info_dict.update({'Nernst-Einstein electrical conductivity': get_std_out([econ, econ_stderr])})

        # diff_gk not use, too slow
        if diff_gk:
            # calculate diffusion constant using Green-Kubo relation
            from ...analyzer.acf import get_t_property_list, get_block_average
            from ...analyzer.fitting import ExpConstfit, ExpConstval
            # fit the data using exponential function
            t_list, diff_list = get_t_property_list(property='diffusion constant', name='System')
            n_block = len([t for t in t_list if t < 1])
            bounds = ([0, 0, 0], [100, 100, 100])
            coef, score = ExpConstfit(get_block_average(t_list, n_block=n_block)[10:],
                                      get_block_average(diff_list * 1e8, n_block=n_block)[10:], bounds=bounds)
            coef[0] /= 1e8
            coef[1] /= 1e8
            diff_gk_dict = {'System': get_std_out([coef[1], ExpConstval(t_list[-1], coef)])}
            for i in range(len(n_mol_list)):
                mol_name = 'MO%i' % (i)
                t_list, diff_list = get_t_property_list(property='diffusion constant', name=mol_name)
                coef, score = ExpConstfit(get_block_average(t_list, n_block=n_block)[10:],
                                          get_block_average(diff_list * 1e8, n_block=n_block)[10:], bounds=bounds)
                coef[0] /= 1e8
                coef[1] /= 1e8
                diff_gk_dict.update({mol_name: get_std_out([coef[1], ExpConstval(t_list[-1], coef)])})
            info_dict.update({'diffusion constant-gk': diff_gk_dict})

            # estimate electrical conductivity using Nernst-Einstein relation
            if charge_list != None and set(charge_list) != {0}:
                econ1 = 0.
                econ2 = 0.
                for i, charge in enumerate(charge_list):
                    mol_name = 'MO%i' % (i)
                    diff1, diff2 = diff_gk_dict.get(mol_name)
                    econ1 += diff1 * charge_list[i] ** 2 * n_mol_list[i]
                    econ2 += diff2 * charge_list[i] ** 2 * n_mol_list[i]
                econ1 *= 1.6 ** 2 / 1.38 * 10 ** 8 / temperature_and_stderr[0] / volume_and_stderr[0]
                econ2 *= 1.6 ** 2 / 1.38 * 10 ** 8 / temperature_and_stderr[0] / volume_and_stderr[0]
                info_dict.update({'Nernst-Einstein electrical conductivity-gk': get_std_out([econ1, econ2])})
            os.remove('traj.gro')

        return info_dict

    # analyze electrical conductivity
    def analyze_econ(self, mstools_dir, weight=0.00):
        from ...panedr import edr_to_df
        df = edr_to_df('nvt.edr')
        temperature = df.Temperature.mean()
        volume = self.gmx.get_volume_from_gro('nvt.gro')
        commands = []
        out, err = self.gmx.current('nvt.trr', 'nvt.tpr', caf=True)
        open('current.out', 'w').write(out)
        open('current.err', 'w').write(err)
        commands.append(os.path.join(mstools_dir, 'mstools', 'cpp', 'current-gk') + ' current.xvg' + ' %f' % (
            volume) + ' %f' % (
                            temperature) + ' %.2f' % (weight))
        for cmd in commands:
            sp = Popen(cmd.split(), stdout=PIPE, stdin=PIPE, stderr=PIPE)
            sp.communicate()

    # analyze viscosity
    def analyze_vis(self,  mstools_dir, weight=0.00):
        df = edr_to_df('nvt.edr')
        temperature = df.Temperature.mean()
        volume = self.gmx.get_volume_from_gro('nvt.gro')
        self.gmx.energy('nvt.edr', properties=['Pres-XY', 'Pres-XZ', 'Pres-YZ'], out='pressure.xvg')
        commands = []
        commands.append(
            os.path.join(mstools_dir, 'mstools', 'cpp', 'vis-gk') + ' pressure.xvg' + ' %f' % (volume) + ' %f' % (
                temperature) + ' %.2f' % (weight))
        for cmd in commands:
            sp = Popen(cmd.split(), stdout=PIPE, stdin=PIPE, stderr=PIPE)
            sp.communicate()

    def analyze_acf(self, mstools_dir, charge_list, n_mol_list, current=False, delete_trr=True, diff_gk=False, weight=0.00):
        info_dict = self.analyze_diff(charge_list, n_mol_list, diff_gk=diff_gk)
        if current:
            # self.analyze_vis(mstools_dir=mstools_dir, weight=weight) this function is implemented in prepare
            self.analyze_econ(mstools_dir=mstools_dir, weight=weight)
        if delete_trr:
            if os.path.exists('nvt.trr'):
                os.remove('nvt.trr')
            if os.path.exists('traj.gro'):
                os.remove('traj.gro')
        info_dict.update({
            'failed': [False],
            'continue': [False],
            'continue_n': 0,
        })

        return info_dict

    @staticmethod
    def post_process(T_list, P_list, result_list, **kwargs) -> (dict, str):
        def round5(x):
            return float('%.5e' % x)
        post_result = {}
        t_set = set(T_list)
        p_set = set(P_list)

        if len(t_set) < 4:
            return None, 'T points less than 4'

        if len(p_set) == 1:
            t_p_viscosity_score_list = []
            t_p_econ_score_list = []
            t_p_NEecon_stderr_list = []
            t_p_diff_list = []
            for i, result in enumerate(result_list):
                if result.get('viscosity') is not None:
                    t_p_viscosity_score_list.append([T_list[i], P_list[i], result.get('viscosity'), result.get('vis_score')]) # [t, p, value, score]
                if result.get('electrical conductivity') is not None:
                    t_p_econ_score_list.append([T_list[i], P_list[i], result.get('electrical conductivity'), result.get('econ_score')]) # [t, p, value, score]
                t_p_NEecon_stderr_list.append([T_list[i], P_list[i], result.get('Nernst-Einstein electrical conductivity')]) # [t, p, [value, stderr]]
                t_p_diff_list.append([T_list[i], P_list[i], result.get('diffusion constant')]) # [t, p, diff_dict{name: [diff, stderr]}]

            t_p_viscosity_score_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_econ_score_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_NEecon_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_diff_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T

            from ...analyzer.fitting import polyfit, VTFfit
            p = str(P_list[0])
            # viscosity VTF fit
            t_vis_VTF = None
            _vis_t_list = [element[0] for element in t_p_viscosity_score_list]
            _vis_list = [element[2] for element in t_p_viscosity_score_list]
            if len(_vis_list) > 6:
                _t_vis_coeff, _t_vis_score = VTFfit(_vis_t_list, _vis_list)
                t_vis_VTF = {}
                t_vis_VTF[p] = [list(map(round5, _t_vis_coeff)), round5(_t_vis_score), min(_vis_t_list),
                                max(_vis_t_list)]
            # electrical conductivity 3th-polyfit
            t_econ_poly3 = None
            _econ_t_list = [element[0] for element in t_p_econ_score_list]
            _econ_list = [element[2] for element in t_p_econ_score_list]
            if len(_econ_list) > 6:
                _t_econ_coeff, _t_econ_score = polyfit(_econ_t_list, _econ_list, 3)
                t_econ_poly3 = {}
                t_econ_poly3[p] = [list(map(round5, _t_econ_coeff)), round5(_t_econ_score), min(_econ_t_list),
                                   max(_econ_t_list)]
            # diffusion constant and Nernst-Einstein electrical conductivity 3th-polyfit
            t_NEecon_poly3 = None
            t_diff_poly3 = None
            _t_list = [element[0] for element in t_p_NEecon_stderr_list]
            _NEecon_list = [element[2][0] for element in t_p_NEecon_stderr_list]
            _name_list = t_p_diff_list[0][2].keys()
            _diff_list = {name: [] for name in _name_list}
            for element in t_p_diff_list:
                for name in _name_list:
                    _diff_list.get(name).append(element[2].get(name)[0])
            if len(_t_list) > 6:
                _t_NEecon_coeff, _t_NEecon_score = polyfit(_t_list, _NEecon_list, 3)
                _t_diff_coeff_score = {}
                for name in _name_list:
                    _t_diff_coeff_score[name] = polyfit(_t_list, _diff_list.get(name), 3)
                t_NEecon_poly3 = {}
                t_diff_poly3 = {}
                t_NEecon_poly3[p] = [list(map(round5, _t_NEecon_coeff)), round5(_t_NEecon_score), min(_t_list),
                                     max(_t_list)]
                t_diff_poly3[p] = [{}, min(_t_list), max(_t_list)]
                for name in _name_list:
                    t_diff_poly3[p][0][name] = [list(map(round5, _t_diff_coeff_score[name][0])),
                                                round5(_t_diff_coeff_score[name][1])]
            post_result = {
                'p': P_list[0],
                'viscosity': t_p_viscosity_score_list, # [t, p, viscosity, score]
                'electrical conductivity': t_p_econ_score_list, # [t, p, econ, score]
                'diffusion constant': t_p_diff_list, # [t, p, diff_dict{name: [diff, stderr]}]
                'Nernst-Einstein electrical conductivity': t_p_NEecon_stderr_list, # [t, p, [value, stderr]]
                'vis-t-VTF': t_vis_VTF, # {'pressure': [[coeff], score, t_min, t_max]}
                'econ-t-poly3': t_econ_poly3,
                'NEecon-t-poly3': t_NEecon_poly3,
                'diff-t-poly3': t_diff_poly3, # {'pressure': [{'System': [[coeff], score]}, t_min, t_max]}
            }
            if result_list[0].get('diffusion constant-gk and score') is not None:
                t_p_diffgk_list = []
                for i, result in enumerate(result_list):
                    t_p_diffgk_list.append([T_list[i], P_list[i], result.get('diffusion constant-gk and score')])
                post_result['diffusion constant gk'] = t_p_diffgk_list # [t, p, diff_dict{name: [diff, score]}]
                '''
                _name_list = t_p_diffgk_list[0][2].keys()
                _diffgk_list = {name: [] for name in _name_list}
                for element in t_p_diffgk_list:
                    for name in _name_list:
                        _diffgk_list.get(name).append(element[2].get(name))
                '''
            return post_result, 'time decomposition method, green-kubo'

        if len(p_set) < 5:
            return None, 'P points less than 5'

    @staticmethod
    def get_post_data(post_result, T, P, **kwargs) -> dict:
        from ...analyzer.fitting import VTFval, polyval

        # single-pressure simulation
        if post_result.get('p') is not None:
            if P != post_result.get('p'):
                raise Exception('for single pressure post_result, P must be the same')

            result = {}
            converge_criterion = 0.95  # R value of fitting
            # viscosity
            if post_result['vis-t-VTF'] is None:
                for t, p, viscosity, score in post_result['viscosity']:
                    if t == T:
                        result['viscosity'] = viscosity
                        break
            else:
                coef, score, tmin, tmax = post_result['vis-t-VTF'][str(P)]
                if score > converge_criterion and T > tmin - 10 and T < tmax + 10:
                    result['viscosity'] = VTFval(T, coef)
            # electrical conductivity
            if post_result['econ-t-poly3'] is None:
                for t, p, electrical_conductivity, score in post_result['electrical conductivity']:
                    if t == T:
                        result['electrical conductivity'] = electrical_conductivity
                        break
            else:
                coef, score, tmin, tmax = post_result['econ-t-poly3'][str(P)]
                if score > converge_criterion and T > tmin - 10 and T < tmax + 10:
                    result['electrical conductivity'] = polyval(T, coef)
            # Nernst-Einstein electrical conductivity
            if post_result['NEecon-t-poly3'] is None:
                for t, p, [electrical_conductivity, stderr] in post_result['Nernst-Einstein electrical conductivity']:
                    if t == T:
                        result['Nernst-Einstein electrical conductivity'] = electrical_conductivity
                        break
            else:
                coef, score, tmin, tmax = post_result['NEecon-t-poly3'][str(P)]
                if score > converge_criterion and T > tmin - 10 and T < tmax + 10:
                    result['Nernst-Einstein electrical conductivity'] = polyval(T, coef)
            # diffusion constant
            if post_result['diff-t-poly3'] is None:
                for t, p, diff in post_result['diffusion constant']:
                    if t == T:
                        result['diffusion constant'] = diff['System'][0]
                        keys = list(diff.keys())
                        keys.remove('System')
                        diff_sum = 0.
                        for i in keys:
                            diff_sum += diff[i][0]
                        if diff_sum != 0.:
                            result['diffusion constant sum'] = diff_sum
                        break
            else:
                diff, tmin, tmax = post_result['diff-t-poly3'][str(P)]
                coef, score = diff['System']
                if score > converge_criterion and T > tmin - 10 and T < tmax + 10:
                    result['diffusion constant'] = polyval(T, coef)
                keys = list(diff.keys())
                keys.remove('System')
                diff_sum = 0.
                for i in keys:
                    coef, score = diff[i]
                    if score > converge_criterion and T > tmin - 10 and T < tmax + 10:
                        diff_sum += polyval(T, coef)
                    else:
                        break
                else:
                    if diff_sum != 0.:
                        result['diffusion constant sum'] = diff_sum
        # multi-pressure simulation
        else:
            # multi-pressure part, need to be finished
            return {}
        return result
