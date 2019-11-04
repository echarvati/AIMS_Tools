import os, shutil
from .gmx import GmxSimulation
from ...wrapper.gmx import *
from ...analyzer.acf import get_std_out

class Nvt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'nvt'
        self.logs = []

    def build(self, ppf=None, minimize=False):
        pass

    def prepare(self, prior_job_dir=None, gro='npt.gro', top='topol.top', T=298, jobname=None,
                dt=0.001, nst_eq=int(4E5), nst_run=int(5E5), random_seed=-1, nst_edr=5, nst_trr=5,
                tcoupl='v-rescale', acf=False, mstools_dir=None, **kwargs) -> [str]:
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

        if acf:
            # diffusion constant
            commands.append(self.gmx.trjconv('nvt.tpr', 'nvt.trr', 'traj.gro', skip=10, get_cmd=True))
            commands.append(os.path.join(mstools_dir, 'mstools', 'cpp', 'diff-gk') + ' traj.gro')
            # viscosity
            commands.append(self.gmx.energy('nvt.edr', properties=['Pres-XY', 'Pres-XZ', 'Pres-YZ'], out='pressure.xvg', get_cmd=True))
            volume = self.gmx.get_volume_from_gro('npt.gro')
            weight = 0.00
            commands.append(os.path.join(mstools_dir, 'mstools', 'cpp', 'vis-gk') + ' pressure.xvg' + ' %f' % (volume) + ' %f' % (T) + ' %.2f' % (weight))
        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands


    # analyze diffusion constant
    def analyze_diff(self, charge_list, n_mol_list):
        # get temperature and volume
        volume_and_stderr = [self.gmx.get_volume_from_gro('nvt.gro'), 0.]
        [temperature_and_stderr] = self.gmx.get_properties_stderr('nvt.edr', ['Temperature'])

        # calculate diffusion constant using Einstein relation
        diff_e_dict = {'System': get_std_out(list(self.gmx.diffusion('nvt.xtc', 'nvt.tpr')))}
        for i in range(len(n_mol_list)):
            mol_name = 'MO%i' % (i)
            diff_e_dict.update({mol_name: get_std_out(list(self.gmx.diffusion('nvt.xtc', 'nvt.tpr', group=mol_name)))})

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
            info_dict.update({'Nernst-Einstein electrical conductivity and standard error via Einstein diffusion constant': get_std_out([econ, econ_stderr])})

        # calculate diffusion constant using Green-Kubo relation
        # os.remove('traj.gro')
        # os.remove('nvt.trr')
        from ...analyzer.acf import get_t_property_list, get_block_average
        from ...analyzer.fitting import ExpConstfit, ExpConstval
        # fit the data using exponential function
        t_list, diff_list = get_t_property_list(property='diffusion constant', name='System')
        n_block = len([t for t in t_list if t < 1])
        coef, score = ExpConstfit(get_block_average(t_list, n_block=n_block)[2:], get_block_average(diff_list, n_block=n_block)[2:])
        diff_gk_dict = {'System': get_std_out([coef[1], ExpConstval(t_list[-1], coef)])}
        for i in range(len(n_mol_list)):
            mol_name = 'MO%i' % (i)
            t_list, diff_list = get_t_property_list(property='diffusion constant', name=mol_name)
            n_block = len([t for t in t_list if t < 1])
            coef, score = ExpConstfit(get_block_average(t_list, n_block=n_block)[2:],
                                      get_block_average(diff_list, n_block=n_block)[2:])
            diff_gk_dict.update({mol_name: get_std_out([coef[1], ExpConstval(t_list[-1], coef)])})
        info_dict.update({'diffusion constant via Green-Kubo relation': diff_gk_dict})

        # estimate electrical conductivity using Nernst-Einstein relation
        if charge_list != None and set(charge_list) != {0}:
            econ1 = 0.
            econ2 = 0.
            for i, charge in enumerate(charge_list):
                mol_name = 'MO%i' % (i)
                diff1, diff2 = diff_gk_dict.get(mol_name)
                econ1 += diff1 * charge_list[i]**2 * n_mol_list[i]
                econ2 += diff2 * charge_list[i]**2 * n_mol_list[i]
            econ1 *= 1.6 ** 2 / 1.38 * 10 ** 8 / temperature_and_stderr[0] / volume_and_stderr[0]
            econ2 *= 1.6 ** 2 / 1.38 * 10 ** 8 / temperature_and_stderr[0] / volume_and_stderr[0]
            info_dict.update({'Nernst-Einstein electrical conductivity and standard error via Green-Kubo diffusion constant': get_std_out([econ1, econ2])})

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

    def analyze_acf(self, mstools_dir, charge_list, n_mol_list, current=False, weight=0.00):
        info_dict = self.analyze_diff(charge_list, n_mol_list)
        if current:
            self.analyze_econ(mstools_dir=mstools_dir, weight=weight)

        info_dict.update({
            'failed': [False],
            'continue': [False],
            'continue_n': 0,
        })

        return info_dict