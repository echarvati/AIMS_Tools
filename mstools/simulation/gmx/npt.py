import os
import shutil, numpy as np
from subprocess import Popen, PIPE

from .gmx import GmxSimulation
from ...analyzer import is_converged, block_average
from ...wrapper.ppf import delta_ppf
from ...panedr import edr_to_df

class Npt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'npt'
        self.logs = ['npt.log', 'hvap.log']
        self.n_atom_default = 3000
        self.n_mol_default = 75
        self.dt = 0.002

    def build(self, export=True, ppf=None):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, self.pdb, size=[i - 2 for i in self.box], silent=True)
        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr=self.pdb, size=self.box)

        # build msd for fast export
        self.packmol.build_box(self.pdb_list, [1] * len(self.pdb_list), self._single_pdb, size=self.box,
                               inp_file='build_single.inp', silent=True)
        self.dff.build_box_after_packmol(self.mol2_list, [1] * len(self.pdb_list), self._single_msd,
                                         mol_corr=self._single_pdb, size=self.box)

        if export:
            self.fast_export_single(ppf=ppf, gro_out='_single.gro', top_out='topol.top')
            self.gmx.pdb2gro(self.pdb, 'conf.gro', [i / 10 for i in self.box], silent=True)  # A to nm
            self.gmx.modify_top_mol_numbers('topol.top', self.n_mol_list)
            if ppf is not None:
                shutil.copy(os.path.join(ppf), 'ff.ppf')

    def prepare(self, model_dir='.', gro='conf.gro', top='topol.top', T=298, P=1, jobname=None, TANNEAL=800,
                dt=0.002, nst_eq=int(4E5), nst_run=int(5E5), nst_edr=100, nst_trr=int(5E4), nst_xtc=int(1E3),
                random_seed=-1, drde=False, tcoupl='langevin', diff_gk=False, mstools_dir=None, **kwargs) -> [str]:
        if os.path.abspath(model_dir) != os.getcwd():
            shutil.copy(os.path.join(model_dir, gro), gro)
            shutil.copy(os.path.join(model_dir, top), top)
            for f in os.listdir(model_dir):
                if f.endswith('.itp'):
                    shutil.copy(os.path.join(model_dir, f), '.')

        if drde:
            ### Temperature dependent parameters
            # TODO Assumes ppf file named ff.ppf
            if os.path.abspath(model_dir) != os.getcwd():
                shutil.copy(os.path.join(model_dir, self._single_msd), self._single_msd)
            delta_ppf(os.path.join(model_dir, 'ff.ppf'), 'ff.ppf', T)
            mol_numbers = self.gmx.get_top_mol_numbers(top)
            self.fast_export_single(ppf='ff.ppf', gro_out='_single.gro', top_out=top)
            self.gmx.modify_top_mol_numbers(top, [n for m, n in mol_numbers])

        nprocs = self.jobmanager.nprocs
        commands = []
        # energy minimization
        self.gmx.prepare_mdp_from_template('t_em.mdp', mdp_out='grompp-em.mdp')
        cmd = self.gmx.grompp(mdp='grompp-em.mdp', gro=gro, top=top, tpr_out='em.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='em', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        gro_em = 'em.gro'
        # NVT annealing from 0 to TANNEAL to target T with Langevin thermostat
        if TANNEAL is not None:
            self.gmx.prepare_mdp_from_template('t_nvt_anneal.mdp', mdp_out='grompp-anneal.mdp', T=T, TANNEAL=TANNEAL,
                                               nsteps=int(1E5), nstxtcout=0)
            cmd = self.gmx.grompp(mdp='grompp-anneal.mdp', gro='em.gro', top=top, tpr_out='anneal.tpr', get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name='anneal', nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

            gro_em = 'anneal.gro'

        # NPT equilibrium with Langevin thermostat and Berendsen barostat
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-eq.mdp', T=T, P=P, gen_seed=random_seed,
                                           nsteps=nst_eq, nstxtcout=0, pcoupl='berendsen')
        cmd = self.gmx.grompp(mdp='grompp-eq.mdp', gro=gro_em, top=top, tpr_out='eq.tpr', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='eq', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # NPT production with Langevin thermostat and Parrinello-Rahman barostat
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-npt.mdp', T=T, P=P,
                                           dt=dt, nsteps=nst_run, nstenergy=nst_edr, nstxout=nst_trr, nstvout=nst_trr,
                                           nstxtcout=nst_xtc, restart=True, tcoupl=tcoupl)
        cmd = self.gmx.grompp(mdp='grompp-npt.mdp', gro='eq.gro', top=top, tpr_out='npt.tpr',
                              cpt='eq.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='npt', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        # Rerun enthalpy of vaporization
        commands.append('export GMX_MAXCONSTRWARN=-1')

        top_hvap = 'topol-hvap.top'
        self.gmx.generate_top_for_hvap(top, top_hvap)
        self.gmx.prepare_mdp_from_template('t_npt.mdp', mdp_out='grompp-hvap.mdp', nstxtcout=0, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-hvap.mdp', gro='eq.gro', top=top_hvap, tpr_out='hvap.tpr', get_cmd=True)
        commands.append(cmd)
        # Use OpenMP instead of MPI when rerun hvap
        cmd = self.gmx.mdrun(name='hvap', nprocs=nprocs, n_omp=nprocs, rerun='npt.xtc', get_cmd=True)
        commands.append(cmd)

        if diff_gk:
            # diffusion constant, do not used, very slow
            commands.append(self.gmx.trjconv('npt.tpr', 'npt.trr', 'traj.gro', skip=10, get_cmd=True))
            commands.append(os.path.join(mstools_dir, 'mstools', 'cpp', 'diff-gk') + ' traj.gro')
        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def extend(self, jobname=None, sh=None, info=None, dt=0.002) -> [str]:
        '''
        if info == none, extend simulation for 500 ps
        '''

        if info is None:
            continue_n = 500 / dt
        elif len(info.get('continue_n')) == 1:
            continue_n = info.get('continue_n')[0]
        else:
            raise Exception('npt.extend(), info.get(\'continue_n\') must be 1 dimensional')
        commands = self.extend_single(jobname=jobname, sh=sh, name=info.get('name')[0], continue_n=continue_n, dt=dt)
        return commands

    def extend_single(self, jobname=None, sh=None, name=None, continue_n=None, dt=0.002) -> [str]:
        '''
        extend simulation for 500 ps
        '''
        nprocs = self.jobmanager.nprocs
        commands = []

        if continue_n is None:
            raise Exception('npt.extend_single(): continue_n cannot be None\n')
        if name != 'npt':
            raise Exception('npt.extend_single(): name must be npt\n')

        extend = continue_n * dt
        self.gmx.extend_tpr('npt.tpr', extend, silent=True)
        # Extending NPT production
        cmd = self.gmx.mdrun(name='npt', nprocs=nprocs, extend=True, get_cmd=True)
        commands.append(cmd)

        # Rerun enthalpy of vaporization
        commands.append('export GMX_MAXCONSTRWARN=-1')
        # Use OpenMP instead of MPI when rerun hvap
        cmd = self.gmx.mdrun(name='hvap', nprocs=nprocs, n_omp=nprocs, rerun='npt.xtc', get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure, sh=sh)
        return commands

    # analyze thermodynamic properties
    def analyze(self, check_converge=True, cutoff_time=7777, **kwargs):
        import numpy as np

        info_dict = {
            'failed': [],
            'continue': [],
            'continue_n': [],
            'reason': [],
            'name': ['npt']
        }
        df = edr_to_df('npt.edr')
        potential_series = df.Potential
        density_series = df.Density
        t_real = df.Temperature.mean()
        length = potential_series.index[-1] # unit: ps, cutoff_time = 7777 ps

        df_hvap = edr_to_df('hvap.edr')
        einter_series = (df_hvap.Potential)

        ### Check the ensemble. KS test on the distribution of kinetic energy
        ### TODO This can be optimized
        import physical_validation as pv
        parser = pv.data.GromacsParser(self.gmx.GMX_BIN)
        data = parser.get_simulation_data(mdp='grompp-npt.mdp', top='topol.top', edr='npt.edr')
        p = pv.kinetic_energy.distribution(data, strict=True, verbosity=0)
        # If test does not pass, set the desired temperature to t_real.
        # Because small deviation in temperature exists for Langevin thermostat
        # 3 Kelvin of deviation is permitted
        if p < 0.01 and abs(data.ensemble.temperature - t_real) < 3:
            try:
                data._SimulationData__ensemble._EnsembleData__t = t_real
                p = pv.kinetic_energy.distribution(data, strict=True, verbosity=0)
            except Exception as e:
                print(repr(e))

        if p < 0.01:
            if length > cutoff_time:
                info_dict.update({'warning: ': 'KS test for kinetic energy failed: p<0.01'})
            else:
                info_dict['failed'].append(False)
                info_dict['continue'].append(True)
                info_dict['continue_n'].append(2.5e5)
                info_dict['reason'].append('KS test for kinetic energy failed: p<0.01')
                return info_dict
        elif p < 0.05:
            if length > cutoff_time:
                info_dict.update({'warning: ': 'KS test for kinetic energy failed: 0.01 < p < 0.05'})
            else:
                info_dict['failed'].append(False)
                info_dict['continue'].append(True)
                info_dict['continue_n'].append(2.5e5)
                info_dict['reason'].append('KS test for kinetic energy failed: 0.01 < p < 0.05')
                return info_dict

        ### Check structure freezing using Density
        if density_series.min() / 1000 < 0.1:  # g/mL
            if length > cutoff_time:
                info_dict['failed'].append(True)
                info_dict['reason'].append('vaporize')
                return info_dict
            else:
                info_dict['failed'].append(False)
                info_dict['continue'].append(True)
                info_dict['continue_n'].append(2.5e5)
                info_dict['reason'].append('vaporize')
                return info_dict

        ### Check structure freezing using Diffusion of COM of molecules. Only use last 400 ps data
        diffusion, _ = self.gmx.diffusion('npt.xtc', 'npt.tpr', mol=True, begin=length - 400)
        if diffusion < 1E-8:  # cm^2/s
            if length > cutoff_time:
                info_dict['failed'].append(True)
                info_dict['reason'].append('freeze')
                return info_dict
            else:
                info_dict['failed'].append(False)
                info_dict['continue'].append(True)
                info_dict['continue_n'].append(2.5e5)
                info_dict['reason'].append('freeze')
                return info_dict

        ### Check convergence
        if check_converge:
            # use potential to do a initial determination
            # use at least 4/5 of the data
            _, when_pe = is_converged(potential_series, frac_min=0)
            when_pe = min(when_pe, length * 0.2)
            # use density to do a final determination
            _, when_dens = is_converged(density_series, frac_min=0)
            when = max(when_pe, when_dens)
            if when > length * 0.5:
                if length > cutoff_time:
                    info_dict.update({'warning: ': 'PE and density not converged'})
                else:
                    info_dict['failed'].append(False)
                    info_dict['continue'].append(True)
                    info_dict['continue_n'].append(2.5e5)
                    info_dict['reason'].append('PE and density not converged')
                    return info_dict
        else:
            when = 0


        ### Get expansion and compressibility using fluctuation method
        nblock = 5
        blocksize = (length - when) / nblock
        expan_list = []
        compr_list = []
        for i in range(nblock):
            begin = when + blocksize * i
            end = when + blocksize * (i + 1)
            expan, compr = self.gmx.get_fluct_props('npt.edr', begin=begin, end=end)
            expan_list.append(expan)
            compr_list.append(compr)
        expansion, expan_stderr = np.mean(expan_list), np.std(expan_list, ddof=1) / np.sqrt(nblock)
        compressi, compr_stderr = np.mean(compr_list), np.std(compr_list, ddof=1) / np.sqrt(nblock)
        expan_stderr = float('%.1e' % expan_stderr)  # 2 effective number for stderr
        compr_stderr = float('%.1e' % compr_stderr)  # 2 effective number for stderr

        temperature_and_stderr, pressure_and_stderr, potential_and_stderr, density_and_stderr, volume_and_stderr, ke_and_stderr, te_and_stderr, pv_and_stderr = \
            self.gmx.get_properties_stderr('npt.edr',
                                           ['Temperature', 'Pressure', 'Potential', 'Density', 'Volume', 'Kinetic-En.', 'Total-Energy', 'pV'],
                                           begin=when)
        if info_dict['failed'] == []:
            info_dict['failed'].append(False)
        if info_dict['continue'] == []:
            info_dict['continue'].append(False)
        if info_dict['reason'] == []:
            info_dict['reason'].append('converge')

        le_and_stderr = []
        le_and_stderr.append(te_and_stderr[0] + pv_and_stderr[0])
        le_and_stderr.append(te_and_stderr[1] + pv_and_stderr[1])
        ad_dict = {
            'density': [i / 1000 for i in density_and_stderr],  # g/mL
            'length'            : length,
            'converge'          : when,
            'temperature'       : temperature_and_stderr,  # K
            'pressure'          : pressure_and_stderr,  # bar
            'potential'         : potential_and_stderr,  # kJ/mol
            'kinetic energy'    : ke_and_stderr, # kJ/mol
            'total energy'      : te_and_stderr, # kJ/mol
            'pV'                : pv_and_stderr,  # kJ/mol
            'liquid enthalpy'   : le_and_stderr, # kJ/mol
            'einter'            : list(block_average(einter_series.loc[when:])),  # kJ/mol
            'expansion'         : [expansion, expan_stderr],
            'compress'          : [compressi, compr_stderr],
        }
        info_dict.update(ad_dict)
        return info_dict

    # analyze diffusion constant
    def analyze_diff(self, charge_list, n_mol_list, diff_gk=False):
        # get temperature and volume
        temperature_and_stderr, volume_and_stderr = self.gmx.get_properties_stderr('npt.edr', ['Temperature', 'Volume'])

        # calculate diffusion constant using Einstein relation
        diff_e_dict = {'System': list(self.gmx.diffusion('npt.xtc', 'npt.tpr'))}
        for i in range(len(n_mol_list)):
            mol_name = 'MO%i' % (i)
            diff_e_dict.update({mol_name: list(self.gmx.diffusion('npt.xtc', 'npt.tpr', group=mol_name))})

        info_dict = {'diffusion constant and standard error via Einstein relation': diff_e_dict}

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
            info_dict.update({'Nernst-Einstein electrical conductivity and standard error via Einstein diffusion constant': [econ, econ_stderr]})

        if diff_gk:
            # calculate diffusion constant using Green-Kubo relation
            from ...analyzer.acf import get_t_property_list, get_block_average
            from ...analyzer.fitting import ExpConstfit, ExpConstval
            # fit the data using exponential function
            t_list, diff_list = get_t_property_list(property='diffusion constant', name='System')
            n_block = len([t for t in t_list if t < 1])
            coef, score = ExpConstfit(get_block_average(t_list, n_block=n_block)[2:],
                                      get_block_average(diff_list, n_block=n_block)[2:])
            diff_gk_dict = {'System': [coef[1], ExpConstval(t_list[-1], coef)]}
            for i in range(len(n_mol_list)):
                mol_name = 'MO%i' % (i)
                t_list, diff_list = get_t_property_list(property='diffusion constant', name=mol_name)
                n_block = len([t for t in t_list if t < 1])
                coef, score = ExpConstfit(get_block_average(t_list, n_block=n_block)[2:],
                                          get_block_average(diff_list, n_block=n_block)[2:])
                diff_gk_dict.update({mol_name: [coef[1], ExpConstval(t_list[-1], coef)]})
            info_dict.update({'diffusion constant via Green-Kubo relation': diff_gk_dict})

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
                info_dict.update({
                                     'Nernst-Einstein electrical conductivity and standard error via Green-Kubo diffusion constant': [
                                         econ1, econ2]})

            os.remove('traj.gro')

        return info_dict

    # analyze electrical conductivity
    def analyze_econ(self, mstools_dir, weight=0.00):
        df = edr_to_df('npt.edr')
        temperature = df.Temperature.mean()
        volume = df.Volume.mean()
        commands = []
        out, err = self.gmx.current('npt.trr', 'npt.tpr', caf=True)
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
        df = edr_to_df('npt.edr')
        temperature = df.Temperature.mean()
        volume = df.Volume.mean()
        self.gmx.energy('npt.edr', properties=['Pres-XY', 'Pres-XZ', 'Pres-YZ'], out='pressure.xvg')
        commands = []
        commands.append(
            os.path.join(mstools_dir, 'mstools', 'cpp', 'vis-gk') + ' pressure.xvg' + ' %f' % (volume) + ' %f' % (
                temperature) + ' %.2f' % (weight))
        for cmd in commands:
            sp = Popen(cmd.split(), stdout=PIPE, stdin=PIPE, stderr=PIPE)
            sp.communicate()

    def analyze_acf(self, mstools_dir, charge_list, n_mol_list, current=False, weight=0.00):
        info_dict = self.analyze_diff(charge_list, n_mol_list)
        self.analyze_vis(mstools_dir=mstools_dir, weight=weight)
        if current:
            self.analyze_econ(mstools_dir=mstools_dir, weight=weight)
        os.remove('npt.trr')
        info_dict.update({
            'failed': [False],
            'continue': [False],
            'continue_n': 0,
        })

        return info_dict

    def clean(self):
        for f in os.listdir(os.getcwd()):
            if f.startswith('em.') or f.startswith('anneal.') or f.startswith('eq.') \
                    or f.endswith('.dfo') or f.startswith('#'):
                try:
                    os.remove(f)
                except:
                    pass

    @staticmethod
    def post_process(T_list, P_list, result_list, n_mol_list, **kwargs) -> (dict, str):
        def round5(x):
            return float('%.5e' % x)
        t_set = set(T_list)
        p_set = set(P_list)
        mols_number = sum(n_mol_list)

        if len(t_set) < 5:
            return None, 'T points less than 5'
        elif len(p_set) == 1:
            dens_stderr_list = [result['density'] for result in result_list]
            eint_stderr_list = [(np.array(result['einter']) / mols_number).tolist() for result in result_list]
            hl_stderr_list = [((np.array(result['einter']) + np.array(result['kinetic energy']) + np.array(result['pV']))
                              / mols_number).tolist() for result in result_list]
            comp_stderr_list = [result['compress'] for result in result_list]

            t_p_dens_stderr_list = list(map(list, zip(T_list, P_list, dens_stderr_list)))
            t_p_dens_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_eint_stderr_list = list(map(list, zip(T_list, P_list, eint_stderr_list)))
            t_p_eint_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_hl_stderr_list = list(map(list, zip(T_list, P_list, hl_stderr_list)))
            t_p_hl_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_comp_stderr_list = list(map(list, zip(T_list, P_list, comp_stderr_list)))
            t_p_comp_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T

            _t_list = [element[0] for element in t_p_dens_stderr_list]
            _dens_list = [element[2][0] for element in t_p_dens_stderr_list]
            _eint_list = [element[2][0] for element in t_p_eint_stderr_list]
            _hl_list = [element[2][0] for element in t_p_hl_stderr_list]
            _comp_list = [element[2][0] for element in t_p_comp_stderr_list]
            from ...analyzer.fitting import polyfit
            _t_dens_coeff, _t_dens_score = polyfit(_t_list, _dens_list, 3)
            _t_eint_coeff, _t_eint_score = polyfit(_t_list, _eint_list, 3)
            _t_hl_coeff, _t_hl_score = polyfit(_t_list, _hl_list, 3)
            _t_comp_coeff, _t_comp_score = polyfit(_t_list, _comp_list, 3)

            p = str(P_list[0])
            t_dens_poly3 = {}
            t_eint_poly3 = {}
            t_hl_poly3 = {}
            t_comp_poly3 = {}
            t_dens_poly3[p] = [list(map(round5, _t_dens_coeff)), round5(_t_dens_score), min(_t_list), max(_t_list)]
            t_eint_poly3[p] = [list(map(round5, _t_eint_coeff)), round5(_t_eint_score), min(_t_list), max(_t_list)]
            t_hl_poly3[p]   = [list(map(round5, _t_hl_coeff  )), round5(_t_hl_score  ), min(_t_list), max(_t_list)]
            t_comp_poly3[p] = [list(map(round5, _t_comp_coeff)), round5(_t_comp_score), min(_t_list), max(_t_list)]

            post_result = {
                'p': P_list[0],
                'density': t_p_dens_stderr_list, # [t, p, [density, stderr]]
                'einter': t_p_eint_stderr_list, # [t, p, [intermolecular energy, stderr]]
                'liquid enthalpy': t_p_hl_stderr_list, # [t, p, [liquid enthalpy, stderr]]
                'compress': t_p_comp_stderr_list, # [t, p, [compressibility, stderr]]
                'dens-t-poly3': t_dens_poly3, # {'pressure': [[coeff], score, t_min, t_max]}
                'einter-t-poly3': t_eint_poly3,
                'hl-t-poly3': t_hl_poly3,
                'compress-t-poly3': t_comp_poly3,
            }
            return post_result, 'Single pressure simulation'
        elif len(p_set) < 5:
            return None, 'P points less than 5'
        else:
            from mstools.analyzer.fitting import polyfit_2d, polyfit, polyval_derivative
            ### einter and liquid enthalpy divided by number of molecules
            dens_stderr_list = [list(map(round5, result['density'])) for result in result_list]
            eint_stderr_list = [list(map(lambda x: round5(x / n_mol_list[0]), result['einter'])) for result in
                                result_list]
            # hl_stderr_list = [list(map(lambda x: round5(x / n_mol_list[0]), result['liquid enthalpy'])) for result in result_list]
            hl_stderr_list = [list(map(lambda x: round5(x / n_mol_list[0]),
                                       [result['einter'][0] + result['kinetic energy'][0] + result['pV'][0],
                                        result['einter'][1] + result['kinetic energy'][1] + result['pV'][1]])) for
                              result in
                              result_list]
            comp_stderr_list = [list(map(round5, result['compress'])) for result in result_list]

            dens_list = [i[0] for i in dens_stderr_list]
            eint_list = [i[0] for i in eint_stderr_list]
            hl_list = [i[0] for i in hl_stderr_list]
            comp_list = [i[0] for i in comp_stderr_list]

            ### Fit with poly4
            coeff_dens, score_dens = polyfit_2d(T_list, P_list, dens_list, 4)
            coeff_eint, score_eint = polyfit_2d(T_list, P_list, eint_list, 4)
            coeff_hl, score_hl = polyfit_2d(T_list, P_list, hl_list, 4)

            ### Fit with poly3
            t_p_dens_list = list(map(list, zip(T_list, P_list, dens_list)))
            t_p_dens_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_eint_list = list(map(list, zip(T_list, P_list, eint_list)))
            t_p_eint_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_hl_list = list(map(list, zip(T_list, P_list, hl_list)))
            t_p_hl_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_comp_list = list(map(list, zip(T_list, P_list, comp_list)))
            t_p_comp_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T

            t_p_dens_stderr_list = list(map(list, zip(T_list, P_list, dens_stderr_list)))
            t_p_dens_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_eint_stderr_list = list(map(list, zip(T_list, P_list, eint_stderr_list)))
            t_p_eint_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_hl_stderr_list = list(map(list, zip(T_list, P_list, hl_stderr_list)))
            t_p_hl_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_comp_stderr_list = list(map(list, zip(T_list, P_list, comp_stderr_list)))
            t_p_comp_stderr_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T

            t_dens_poly3 = {}
            t_eint_poly3 = {}
            t_hl_poly3 = {}
            t_comp_poly3 = {}

            for p in sorted(p_set):
                _t_list = [element[0] for element in t_p_dens_list if element[1] == p]
                _dens_list = [element[2] for element in t_p_dens_list if element[1] == p]
                _eint_list = [element[2] for element in t_p_eint_list if element[1] == p]
                _hl_list = [element[2] for element in t_p_hl_list if element[1] == p]
                _comp_list = [element[2] for element in t_p_comp_list if element[1] == p]

                if len(_t_list) < 5:
                    continue

                # density-T relation is fitted by 3-order polynomial function
                _t_dens_coeff, _t_dens_score = polyfit(_t_list, _dens_list, 3)
                _t_eint_coeff, _t_eint_score = polyfit(_t_list, _eint_list, 3)
                _t_hl_coeff, _t_hl_score = polyfit(_t_list, _hl_list, 3)
                _t_comp_coeff, _t_comp_score = polyfit(_t_list, _comp_list, 3)

                t_dens_poly3[p] = [list(map(round5, _t_dens_coeff)), round5(_t_dens_score), min(_t_list), max(_t_list)]
                t_eint_poly3[p] = [list(map(round5, _t_eint_coeff)), round5(_t_eint_score), min(_t_list), max(_t_list)]
                t_hl_poly3[p] = [list(map(round5, _t_hl_coeff)), round5(_t_hl_score), min(_t_list), max(_t_list)]
                t_comp_poly3[p] = [list(map(round5, _t_comp_coeff)), round5(_t_comp_score), min(_t_list), max(_t_list)]

            post_result = {
                'density': t_p_dens_stderr_list,  # [t, p, [dens]]
                'einter': t_p_eint_stderr_list,
                'liquid enthalpy': t_p_hl_stderr_list,
                'compress': t_p_comp_stderr_list,
                'density-poly4': [list(map(round5, coeff_dens)), round5(score_dens)],
                'einter-poly4': [list(map(round5, coeff_eint)), round5(score_eint)],
                'hl-poly4': [list(map(round5, coeff_hl)), round5(score_hl)],
                'dens-t-poly3': t_dens_poly3,
                'einter-t-poly3': t_eint_poly3,
                'hl-t-poly3': t_hl_poly3,
                'compress-t-poly3': t_comp_poly3,
            }

            return post_result, 'density-poly4-score %.4f einter-poly4-score %.4f' % (score_dens, score_eint)

    @staticmethod
    def get_post_data(post_result, T, P, smiles_list, **kwargs) -> dict:
        from mstools.analyzer.fitting import polyval_derivative_2d, polyval, polyval_derivative, polyfit
        import pybel

        ### Default value
        result = {}
        molwt = 0.
        converge_criterion = 0.95 # R value of fitting
        for smiles in smiles_list:
            py_mol = pybel.readstring('smi', smiles)
            molwt += py_mol.molwt

        # single-pressure simulation
        if post_result.get('p') is not None:
            if P != post_result.get('p'):
                raise Exception('for single pressure post_result, P must be the same')

            # density
            coef, score, tmin, tmax = post_result['dens-t-poly3'][str(P)]
            if score > converge_criterion and T > tmin - 10 and T < tmax + 10:
                density, dDdT = polyval_derivative(T, coef)
                result['density'] = density
                result['expansion'] = -1 / density * dDdT  # K^-1
                result['cp_pv'] = - molwt * P / density ** 2 * dDdT * 0.1  # J/mol/K
            # einter
            coef, score, tmin, tmax = post_result['einter-t-poly3'][str(P)]
            if score > converge_criterion and T > tmin - 10 and T < tmax + 10:
                einter, dEdT = polyval_derivative(T, coef)
                result['einter'] = einter
                result['hvap'] = 8.314 * T / 1000 - einter  # kJ/mol
                result['cp_inter'] = dEdT * 1000  # J/mol.K
            # compressibility
            coef, score, tmin, tmax = post_result['compress-t-poly3'][str(P)]
            if score > converge_criterion and T > tmin - 10 and T < tmax + 10:
                result['compressibility'] = polyval(T, coef)
        # multi-pressure simulation
        else:
            ### Calculate with T,P-poly4. Not accurate enough, especially for expansion and compressibility
            if post_result.get('density-poly4') != None:
                coeff_dens, score_dens = post_result['density-poly4']
                density4, dDdT4, dDdP4 = polyval_derivative_2d(T, P, 4, coeff_dens)  # g/mL
                expansion4 = -1 / density4 * dDdT4  # K^-1
                compressibility4 = 1 / density4 * dDdP4  # bar^-1
                cp_pv4 = - molwt * P / density4 ** 2 * dDdT4 * 0.1  # J/mol/K
                ad_dict = {
                    'density-poly4-score': score_dens,
                    'density-poly4': density4,
                    'expansion-poly4': expansion4,
                    'compress-poly4': compressibility4,
                    'cp_pv-poly4': cp_pv4,
                }
                result.update(ad_dict)
            if post_result.get('einter-poly4') != None:
                coeff_eint, score_eint = post_result['einter-poly4']
                einter4, dEdT4, dEdP4 = polyval_derivative_2d(T, P, 4, coeff_eint)  # kJ/mol

                cp_inter4 = dEdT4 * 1000  # J/mol.K
                ad_dict = {
                    'einter-poly4-score': score_eint,
                    'einter-poly4': einter4,
                    'cp_inter-poly4': cp_inter4,
                }
                result.update(ad_dict)
            if len(post_result['dens-t-poly3']) >= 5:
                _p_dens_list = []
                _p_dDdT_list = []
                for _p in post_result['dens-t-poly3']:
                    coef, score, tmin, tmax = post_result['dens-t-poly3'][str(_p)]
                    if score < converge_criterion or T < tmin - 10 or T > tmax + 10:
                        continue

                    dens, dDdT = polyval_derivative(T, coef)
                    _p_dens_list.append([float(_p), dens])
                    _p_dDdT_list.append([float(_p), dDdT])
                if len(_p_dens_list) >= 5:
                    coef, score = polyfit(*zip(*_p_dens_list), 3)
                    _p_list = list(zip(*_p_dens_list))[0]
                    if P > min(_p_list) - 10 and P < max(_p_list) + 10:
                        density = polyval(P, coef)
                        coef, score = polyfit(*zip(*_p_dDdT_list), 3)
                        dDdT = polyval(P, coef)
                        result['density'] = density
                        result['expansion'] = -1 / density * dDdT  # K^-1
                        result['cp_pv'] = - molwt * P / density ** 2 * dDdT * 0.1  # J/mol/K
            if len(post_result['einter-t-poly3']) >= 5:
                _p_eint_list = []
                _p_dEdT_list = []
                for _p in post_result['einter-t-poly3']:
                    coef, score, tmin, tmax = post_result['einter-t-poly3'][str(_p)]
                    if score < converge_criterion or T < tmin - 10 or T > tmax + 10:
                        continue

                    eint, dEdT = polyval_derivative(T, coef)
                    _p_eint_list.append([float(_p), eint])
                    _p_dEdT_list.append([float(_p), dEdT])
                if len(_p_eint_list) >= 5:
                    coef, score = polyfit(*zip(*_p_eint_list), 3)
                    _p_list = list(zip(*_p_eint_list))[0]
                    if P > min(_p_list) - 10 and P < max(_p_list) + 10:
                        einter = polyval(P, coef)
                        result['einter'] = einter
                        result['hvap'] = 8.314 * T / 1000 - einter  # kJ/mol

                        coef, score = polyfit(*zip(*_p_dEdT_list), 3)
                        dEdT = polyval(P, coef)
                        result['cp_inter'] = dEdT * 1000  # J/mol.K
            if len(post_result['compress-t-poly3']) >= 5:
                _p_comp_list = []
                for _p in post_result['compress-t-poly3']:
                    coef, score, tmin, tmax = post_result['compress-t-poly3'][str(_p)]
                    if score < converge_criterion or T < tmin - 10 or T > tmax + 10:
                        continue

                    _p_comp_list.append([float(_p), polyval(T, coef)])
                if len(_p_comp_list) >= 5:
                    coef, score = polyfit(*zip(*_p_comp_list), 3)
                    _p_list = list(zip(*_p_comp_list))[0]
                    if P > min(_p_list) - 10 and P < max(_p_list) + 10:
                        result['compressibility'] = polyval(P, coef)
            if len(post_result['hl-t-poly3']) >= 5:
                _p_hl_list = []
                for _p in post_result['hl-t-poly3']:
                    coef, score, tmin, tmax = post_result['hl-t-poly3'][str(_p)]
                    if score < converge_criterion or T < tmin - 10 or T > tmax + 10:
                        continue

                    hl = polyval(T, coef)
                    _p_hl_list.append([float(_p), hl])
                if len(_p_hl_list) >= 5:
                    coef, score = polyfit(*zip(*_p_hl_list), 3)
                    _p_list = list(zip(*_p_hl_list))[0]
                    if P > min(_p_list) - 10 and P < max(_p_list) + 10:
                        result['liquid_enthalpy'] = polyval(P, coef)  # kJ/mol

        return result
