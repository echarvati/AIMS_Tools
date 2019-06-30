import os
import shutil
import math
from collections import OrderedDict

from .gmx import GmxSimulation
from ...analyzer import is_converged, block_average, average_of_blocks
from ..trajectory import Trajectory

class NptPPM(GmxSimulation):
    def __init__(self, amplitudes_steps=None, **kwargs):
        super().__init__(**kwargs)
        self.procedure = 'npt-ppm'
        self.dt = 0.001
        self.n_atoms_default = 6000
        self.amplitudes_steps = amplitudes_steps or OrderedDict([(0.010, int(6.0e6)),
                                                                 (0.020, int(3.0e6)),
                                                                 (0.030, int(2.0e6)),
                                                                 (0.040, int(2.0e6)),
                                                                 # (0.050, int(1.0e6)),
                                                                 ])
        self.logs = ['ppm-%.3f.log' % ppm for ppm in self.amplitudes_steps.keys()]

    def build(self, export=True, ppf=None):
        print('Build coordinates using Packmol: %s molecules ...' % self.n_mol_list)
        self.packmol.build_box(self.pdb_list, self.n_mol_list, 'init.pdb', length=self.length - 2, silent=True)

        print('Create box using DFF ...')
        self.dff.build_box_after_packmol(self.mol2_list, self.n_mol_list, self.msd, mol_corr='init.pdb',
                                         length=self.length)
        if export:
            self.export(ppf=ppf)

    def prepare(self, prior_job_dir=None, gro='conf.gro', top='topol.top', T=298, P=1, jobname=None,
                dt=0.001, nst_eq=int(1E5), nst_edr=50, replicate=None,
                **kwargs) -> [str]:
        self.dt = dt
        if prior_job_dir is None:
            raise Exception('prior_job_dir is needed for PPM simulation')

        # Copy gro and topology files from prior NPT simulation
        shutil.copy(os.path.join(prior_job_dir, gro), '.')
        shutil.copy(os.path.join(prior_job_dir, top), '.')
        for f in os.listdir(prior_job_dir):
            if f.endswith('.itp'):
                shutil.copy(os.path.join(prior_job_dir, f), '.')

        if replicate is not None:
            self.gmx.replicate_gro(gro, top, replicate)

        nprocs = self.jobmanager.nprocs
        commands = []
        '''
        if restart:
            for ppm, nst_run in self.amplitudes_steps.items():
                dt = nst_run * 0.001# ppm simulation use 1fs timestep
                name_ppm = 'ppm-%.3f' % ppm
                self.gmx.extend_tpr(name_ppm + '.tpr', dt, silent=True)
                cmd = self.gmx.mdrun(name=name_ppm, nprocs=nprocs, extend=True, get_cmd=True)
                commands.append(cmd)
        '''
        for ppm, nst_run in self.amplitudes_steps.items():
            name_eq = 'eq-%.3f' % ppm
            name_ppm = 'ppm-%.3f' % ppm

            # NPT-PPM equilibrium with Nose-Hoover thermostat and Parrinello-Rahman barostat
            # TODO should test the validity of V-rescale thermostat. V-rescale is preferred if it works
            self.gmx.prepare_mdp_from_template('t_npt_ppm.mdp', mdp_out='grompp-%s.mdp' % name_eq, T=T, P=P,
                                               nsteps=nst_eq, nstxtcout=0, restart=True,
                                               tcoupl='v-rescale', ppm=ppm)
            cmd = self.gmx.grompp(mdp='grompp-%s.mdp' % name_eq, gro=gro, top=top,
                                  tpr_out='%s.tpr' % name_eq, get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name_eq, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

            # NPT-PPM production with Nose-Hoover thermostat and Parrinello-Rahman barostat
            # TODO should test the validity of V-rescale thermostat. V-rescale is preferred if it works
            self.gmx.prepare_mdp_from_template('t_npt_ppm.mdp', mdp_out='grompp-%s.mdp' % name_ppm, T=T, P=P,
                                               dt=dt, nsteps=nst_run, nstenergy=nst_edr, restart=True,
                                               tcoupl='v-rescale', ppm=ppm, nstxtcout=100000)
            cmd = self.gmx.grompp(mdp='grompp-%s.mdp' % name_ppm, gro='%s.gro' % name_eq, top=top,
                                  cpt='%s.cpt' % name_eq, tpr_out='%s.tpr' % name_ppm, get_cmd=True)
            commands.append(cmd)
            cmd = self.gmx.mdrun(name=name_ppm, nprocs=nprocs, get_cmd=True)
            commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure)
        return commands

    def extend(self, jobname=None, sh=None, info=None, dt=0.001) -> [str]:
        '''
        extend simulation for 500 ps
        '''
        nprocs = self.jobmanager.nprocs
        commands = []
        for i, tpr in enumerate(info.get('name')):
            if info.get('continue')[i]:
                extend = info.get('continue_n')[i] * dt
                self.gmx.extend_tpr('%s.tpr' % tpr, extend, silent=True)
                # Extending PPM production
                cmd = self.gmx.mdrun(name=tpr, nprocs=nprocs, extend=True, get_cmd=True)
                commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure, sh=sh)
        return commands

    def ppm_is_converged(self, trj):
        traj = Trajectory(trj, readmass=True, COM=True)
        frame = traj.traj_info[-1]
        frame.ReducedUnitTransform()
        nst = frame.step - traj.traj_info[0].step
        n1 = n2 = n3 = n = 0
        X1 = Y1 = Z1 = M1 = 0.
        X2 = Y2 = Z2 = M2 = 0.
        X3 = Y3 = Z3 = M3 = 0.
        for i in range(frame.atom_number):
            # print(frame.mass[i])
            if (not 0<frame.x[i]<1) or (not 0<frame.y[i]<1) or (not 0<frame.z[i]<1):
                n += 1
            if frame.x[i] > 1:
                if frame.z[i] > 0.5:
                    X1 += frame.x[i] * frame.mass[i]
                    Y1 += frame.y[i] * frame.mass[i]
                    Z1 += frame.z[i] * frame.mass[i]
                    M1 += frame.mass[i]
                    # print(M1)
                    n1 += 1
                else:
                    X2 += frame.x[i] * frame.mass[i]
                    Y2 += frame.y[i] * frame.mass[i]
                    Z2 += frame.z[i] * frame.mass[i]
                    M2 += frame.mass[i]
                    n2 += 1
            elif frame.x[i] < 0:
                X3 += frame.x[i] * frame.mass[i]
                Y3 += frame.y[i] * frame.mass[i]
                Z3 += frame.z[i] * frame.mass[i]
                M3 += frame.mass[i]
                n3 += 1
        if n1 == 0 or n2 == 0 or n3 == 0:
            if nst > 1.4e7:
                return {
                'failed': True,
                'reason': 'the viscosity of this liquid is too high',
                'continue': False,
                'continue_n': 0
            }
            else:
                return {
                'failed': False,
                'reason': None,
                'continue': True,
                'continue_n': int(9.0e6)
            }
        X1 /= M1; Y1 /= M1; Z1 /= M1; n1 /= frame.atom_number
        X2 /= M2; Y2 /= M2; Z2 /= M2; n2 /= frame.atom_number
        X3 /= M3; Y3 /= M3; Z3 /= M3; n3 /= frame.atom_number
        n /= frame.atom_number
        if n > 0.75:
            return{
                'failed': False,
                'reason': 'converged',
                'continue': False,
                'continue_n': 0
            }
        if X1 > 1.6:
            if 0.8 < Z1 < 1.1 and n1 > 0.15:
                rst1 = 0

            else:
                return {
                    'failed': True,
                    'reason': 'PPM trajectory error, please check it manually'
                }
        else:
            rst1 = nst / (X1-1) * 0.6 - nst
        if X2 > 1.6:
            if -0.1 < Z2 < 0.2 and n2 > 0.15:
                rst2 = 0
            else:
                return {
                    'failed': True,
                    'reason': 'PPM trajectory error, please check it manually'
                }
        else:
            rst2 = nst / (X1-1) * 0.6 - nst
        if X3 < -0.6:
            if 0.45 < Z3 < 0.55 and n3 > 0.3:
                rst3 = 0
            else:
                return {
                    'failed': True,
                    'reason': 'PPM trajectory error, please check it manually'
                }
        else:
            rst3 = nst / (X1-1) * 0.6 - nst
        rst = max(rst1, rst2, rst3)
        rst = math.ceil(rst / 1.0e6) * 1.0e6
        if rst > 1e8 and nst > 1.4e7:
            return {
                'failed': True,
                'reason': 'the viscosity of this liquid is too high, need approximately %i additional step to converge' % (rst),
                'continue': False,
                'continue_n': 0
            }
        result = {
            'failed': False,
            'reason': 'not converged',
            'continue': True,
            'continue_n': min(rst, int(9.0e6))
        }
        if rst == 0:
            result['continue'] = False
            result['reason'] = 'converged'
        return result

    def analyze(self, dirs=None, check_converge=True, prior_result=None, **kwargs):
        import numpy as np
        from ...panedr import edr_to_df
        from ...analyzer.fitting import polyfit

        vis_list = []
        stderr_list = []
        info_dict = {
            'failed': [],
            'continue': [],
            'reason': [],
            'name': [],
            'continue_n': [],
            'warning': []
        }
        for ppm in self.amplitudes_steps.keys():
            name_ppm = 'ppm-%.3f' % ppm
            df = edr_to_df('%s.edr' % name_ppm)

            # density check
            density_series = df.Density
            potential_series = df.Potential
            length = potential_series.index[-1]
            ### Check structure freezing using Density
            if density_series.min() / 1000 < 0.1:  # g/mL
                info_dict['failed'].append(True)
                info_dict['continue'].append(None)
                info_dict['continue_n'].append(None)
                info_dict['reason'].append('vaporize')
                info_dict['name'].append(name_ppm)
            ### Check convergence
            else:
                if check_converge:
                    # use potential to do a initial determination
                    # use at least 4/5 of the data
                    _, when_pe = is_converged(potential_series, frac_min=0)
                    when_pe = min(when_pe, length * 0.2)
                    # use density to do a final determination
                    _, when_dens = is_converged(density_series, frac_min=0)
                    when = max(when_pe, when_dens)
                    if when > length * 0.5:
                        info_dict['failed'].append(False)
                        info_dict['continue'].append(True)
                        info_dict['continue_n'].append(self.amplitudes_steps[ppm])
                        info_dict['reason'].append('PE and density not converged')
                        info_dict['name'].append(name_ppm)
                    else:
                        if prior_result != None:
                            density_and_stderr = self.gmx.get_properties_stderr('%s.edr' % (name_ppm), ['Density'], begin=when)
                            den_dev = abs(density_and_stderr[0][0] / 1000 - prior_result.get('density')[0]) / prior_result.get('density')[0]
                            if den_dev > 0.02:
                                info_dict['warning'].append('density is inconsistent with prior npt simulation %f, %f' % (prior_result.get('density')[0], density_and_stderr[0][0]/1000))
                            else:
                                info_dict['warning'].append(None)

                        self.gmx.trjconv('%s.tpr' % name_ppm, '%s.xtc' % name_ppm, '%s_trj.gro' % name_ppm,
                                         pbc_nojump=True, silent=True)
                        result = self.ppm_is_converged('%s_trj.gro' % name_ppm)
                        info_dict['failed'].append(result.get('failed'))
                        info_dict['continue'].append(result.get('continue'))
                        info_dict['continue_n'].append(result.get('continue_n'))
                        info_dict['reason'].append(result.get('reason'))
                        info_dict['name'].append(name_ppm)
                        os.remove('%s_trj.gro' % name_ppm)

            ###

            inv_series = df['1/Viscosity']

            # select last 4/5 of data
            when = inv_series.index[len(inv_series) // 5]

            # use block average to estimate stderr, because 1/viscosity fluctuate heavily
            inv_blocks = average_of_blocks(inv_series.loc[when:])
            vis_blocks = [1000 / inv for inv in inv_blocks]  # convert Pa*s to cP
            vis_list.append(np.mean(vis_blocks))
            stderr_list.append(np.std(vis_blocks, ddof=1) / math.sqrt(len(vis_blocks)))

        if set(info_dict.get('failed'))=={False} and set(info_dict.get('continue'))=={False}:
            coef_, score = polyfit(self.amplitudes_steps.keys(), vis_list, 1)
        # coef_, score = polyfit(self.amplitudes_steps.keys(), vis_list, 1, weight=1 / np.sqrt(stderr_list))
            ad_dict = {
            'viscosity'  : coef_[0],
            'score'      : score,
            'vis_list'   : vis_list,
            'stderr_list': stderr_list,
            }
            info_dict.update(ad_dict)

        return info_dict

    def clean(self):
        for f in os.listdir(os.getcwd()):
            if f.startswith('eq.') or f.startswith('#'):
                try:
                    os.remove(f)
                except:
                    pass
