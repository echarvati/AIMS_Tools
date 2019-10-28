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
                dt=0.001, nst_eq=int(4E5), nst_run=int(5E5), random_seed=-1, nstenergy_vis=5, nst_trr=int(5E4),
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
                                           nstenergy=nstenergy_vis, nstxout=nst_trr, nstvout=nst_trr, restart=True)
        cmd = self.gmx.grompp(mdp='grompp-nvt.mdp', gro='eq.gro', top=top, tpr_out='nvt.tpr',
                              cpt='eq.cpt', get_cmd=True)
        commands.append(cmd)
        cmd = self.gmx.mdrun(name='nvt', nprocs=nprocs, get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def analyze(self, skip=1, current=False, mstools_dir=None, temperature=None):
        if mstools_dir is None:
            return {
                'failed': [True],
                'continue': [False],
                'continue_n': 0,
                'reason': 'nvt.analyze(mstools_dir=None) mstools_dir cannot be None'
            }
        if temperature is None:
            return {
                'failed': [True],
                'continue': [False],
                'continue_n': 0,
                'reason': 'nvt.analyze(T=None) T cannot be None'
            }

        from subprocess import Popen, PIPE
        self.gmx.energy('nvt.edr', properties=['Pres-XY', 'Pres-XZ', 'Pres-YZ'], skip=skip, out='pressure.xvg')
        volume = self.gmx.get_volume_from_gro('nvt.gro')
        commands = [
            os.path.join(mstools_dir, 'mstools', 'cpp', 'vis-gk') + ' pressure.xvg' + ' %s' % (volume) + ' %s' % (
                temperature)]
        if current:
            self.gmx.current(acf=True)
            commands.append(os.path.join(mstools_dir, 'mstools', 'cpp', 'current-gk') + 'acf.xvg')

        for cmd in commands:
            sp = Popen(cmd.split(), stdout=PIPE, stdin=PIPE, stderr=PIPE)
            sp.communicate()
        info_dict = {
            'failed': [False],
            'continue': [False],
            'continue_n': 0,
        }

        return info_dict
