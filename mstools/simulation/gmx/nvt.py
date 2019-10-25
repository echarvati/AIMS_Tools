import os, shutil
from .gmx import GmxSimulation


class Nvt(GmxSimulation):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'nvt'
        self.logs = []

    def build(self, ppf=None, minimize=False):
        pass

    def prepare(self, prior_job_dir=None, gro='npt.gro', top='topol.top', T=298, jobname=None,
                dt=0.001, nst_eq=int(4E5), nst_run=int(5E5), random_seed=-1, nstenergy_vis=5, nstvout=5,
                tcoupl='v-rescale', **kwargs) -> [str]:
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
                                           nstenergy=nstenergy_vis, nstvout=nstvout, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-nvt.mdp', gro='eq.gro', top=top, tpr_out='nvt.tpr',
                              cpt='eq.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='nvt', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def analyze(self, skip=1):
        from ...analyzer.acf import get_acf, get_integral
        from ...panedr import edr_to_df
        def get_skip_data(datas, skip=1):
            result = []
            for i, data in enumerate(datas):
                if i % skip == 0:
                    result.append(data)
            return result
        df = edr_to_df('nvt.edr')
        time = df['Time'].tolist()
        pxy = get_skip_data(df['Pres-XY'], skip=skip)
        pxz = get_skip_data(df['Pres-XZ'], skip=skip)
        pyz = get_skip_data(df['Pres-YZ'], skip=skip)
        t_real = df.Temperature.mean()
        V_real = self.gmx.get_box_from_gro('nvt.gro')
        a1, b1 = get_acf(time, pxy)
        a2, b2 = get_acf(time, pxz)
        a3, b3 = get_acf(time, pyz)
        a, b = get_integral(a1, (b1 + b2 + b3) / 3)
        info_dict = {
            'failed': [False],
            'continue': [False],
            'continue_n': 0,
            't_list': a.tolist(),
            'vis_list': (b * 6.022 *10**(-3) * V_real / (8.314 * t_real)).tolist(),
        }
        return info_dict
