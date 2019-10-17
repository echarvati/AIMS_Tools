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
        self.amplitudes_steps = amplitudes_steps or OrderedDict([(0.010, int(4.5e6)),
                                                                 (0.020, int(3.0e6)),
                                                                 (0.030, int(2.0e6)),
                                                                 (0.040, int(1.0e6)),
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
                dt=0.001, nst_eq=int(1E5), nst_edr=500, replicate=None,
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

    # extend function used when EXTEND_GMX_MULTI = False
    def extend(self, jobname=None, sh=None, info=None, dt=0.001) -> [str]:
        '''
        extend simulation for 500 ps
        '''
        if info == None:
            raise Exception('npt_ppm.extend(): info cannot be none')
        nprocs = self.jobmanager.nprocs
        commands = []
        for i, name in enumerate(info.get('name')):
            if info.get('continue')[i]:
                extend = info.get('continue_n')[i] * dt
                self.gmx.extend_tpr('%s.tpr' % name, extend, silent=True)
                # Extending PPM production
                cmd = self.gmx.mdrun(name=name, nprocs=nprocs, extend=True, get_cmd=True)
                commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure, sh=sh)
        return commands

    # extend function used when EXTEND_GMX_MULTI = True
    # def extend_multi(self):
    def extend_single(self, jobname=None, sh=None, name=None, continue_n=None, dt=0.001) -> [str]:
        '''
        extend simulation for 500 ps
        '''
        if continue_n == None:
            raise Exception('npt_ppm.extend_single(): continue_n cannot be None\n')
        if name == None:
            raise Exception('npt_ppm.extend_single(): name cannot be None\n')
        nprocs = self.jobmanager.nprocs
        commands = []
        extend = continue_n * dt
        self.gmx.extend_tpr('%s.tpr' % name, extend, silent=True)
        # Extending PPM production
        cmd = self.gmx.mdrun(name=name, nprocs=nprocs, extend=True, get_cmd=True)
        commands.append(cmd)

        self.jobmanager.generate_sh(os.getcwd(), commands, name=jobname or self.procedure, sh=sh)
        return commands

    def ppm_is_converged(self, trj):
        traj = Trajectory(trj, readmass=True, COM=True, head_and_tail=True)
        frame = traj.traj_info[-1]
        frame.ReducedUnitTransform()
        nst = frame.step - traj.traj_info[0].step
        if nst == 0:
            nst = (frame.t - traj.traj_info[0].t) / self.dt
        if nst == 0:
            return {
                'failed': True,
                'reason': 'The simulation ended abnormally in this case, mostly occur in long chain system due to larger acceleration',
                'continue': False,
                'continue_n': 0
            }
        n1 = n2 = n3 = n = 0
        X1 = Y1 = Z1 = M1 = 0.
        X2 = Y2 = Z2 = M2 = 0.
        X3 = Y3 = Z3 = M3 = 0.
        for i in range(frame.atom_number):
            if (not 0<frame.x[i]<1) or (not 0<frame.y[i]<1) or (not 0<frame.z[i]<1):
                n += 1
            if frame.x[i] > 1:
                if frame.z[i] > 0.5:
                    X1 += frame.x[i] * frame.mass[i]
                    Y1 += frame.y[i] * frame.mass[i]
                    Z1 += frame.z[i] * frame.mass[i]
                    M1 += frame.mass[i]
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
            if nst > 2.0e7:
                return {
                'failed': True,
                'reason': 'the viscosity of this liquid is too high, set this job failed',
                'continue': False,
                'continue_n': 0
            }
            else:
                return {
                'failed': False,
                'reason': 'the viscosity of this liquid is too high, continue this job',
                'continue': True,
                'continue_n': int(1.0e7)
            }
        X1 /= M1; Y1 /= M1; Z1 /= M1; n1 /= frame.atom_number
        X2 /= M2; Y2 /= M2; Z2 /= M2; n2 /= frame.atom_number
        X3 /= M3; Y3 /= M3; Z3 /= M3; n3 /= frame.atom_number
        n /= frame.atom_number
        ad_dict = {
            'more_info': 'X1 = %f, Y1 = %f, Z1 = %f, n1 = %f, X2 = %f, Y2 = %f, Z2 = %f, n2 = %f, X3 = %f, Y3 = %f, Z3 = %f, n3 = %f' % (X1, Y1, Z1, n1, X2, Y2, Z2, n2, X3, Y3, Z3, n3)
        }
        if n > 0.75:
            return{
                'failed': False,
                'reason': 'converged',
                'continue': False,
                'continue_n': 0
            }
        converge_criterion = 0.3
        if X1 > 1 + converge_criterion:
            rst1 = 0
        elif X1 < 1:
            raise Exception('X1 should > 1, unknown error occur')
        else:
            rst1 = nst / (X1-1) * converge_criterion - nst
        if X2 > 1 + converge_criterion:
            rst2 = 0
        elif X2 < 1:
            raise Exception('X2 should > 1, unknown error occur')
        else:
            rst2 = nst / (X2-1) * converge_criterion - nst
        if X3 < -converge_criterion:
            rst3 = 0
        elif X3 > 0:
            raise Exception('X3 should < 0, unknown error occur')
        else:
            rst3 = nst / (-X3) * converge_criterion - nst
        rst = max(rst1, rst2, rst3)
        rst = int(math.ceil(rst / 1.0e6) * 1.0e6)
        if rst > 5e8 and nst > 2.0e7:
            info_dict = {
                'failed': True,
                'reason': 'the viscosity of this liquid is too high, need approximately %i additional step to converge' % (rst),
                'continue': False,
                'continue_n': 0
            }
        else:
            info_dict = {
                'failed': False,
                'reason': 'not converged',
                'continue': True,
                'continue_n': min(rst, int(1.0e7))
            }
        if rst == 0:
            info_dict['continue'] = False
            info_dict['reason'] = 'converged'
        info_dict.update(ad_dict)
        return info_dict

    def analyze(self, dirs=None, check_converge=True, more_info=True, **kwargs):
        import numpy as np
        from ...panedr import edr_to_df
        from ...analyzer.fitting import polyfit

        a_list = []
        vis_list = []
        stderr_list = []
        info_dict = {
            'failed': [],
            'continue': [],
            'continue_n': [],
            'reason': [],
            'name': [],
            'warning': [],
            'more_info': [],
            'length': []
        }
        for ppm in self.amplitudes_steps.keys():
            name_ppm = 'ppm-%.3f' % ppm
            if not os.path.exists('%s.edr' % name_ppm):
                info_dict['failed'].append(True)
                info_dict['continue'].append(False)
                info_dict['continue_n'].append(0)
                info_dict['reason'].append('file do not exists')
                continue

            df = edr_to_df('%s.edr' % name_ppm)

            # density check
            density_series = df.Density
            potential_series = df.Potential
            length = potential_series.index[-1]

            info_dict['name'].append(name_ppm)
            info_dict['length'].append(length)
            ### Check structure freezing using Density
            if density_series.min() / 1000 < 0.1:  # g/mL
                info_dict['failed'].append(True)
                info_dict['continue'].append(False)
                info_dict['continue_n'].append(0)
                info_dict['reason'].append('vaporize')
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
                        info_dict['warning'].append(None)
                        if more_info:
                            info_dict['more_info'].append(None)
                    else:
                        self.gmx.trjconv('%s.tpr' % name_ppm, '%s.xtc' % name_ppm, '%s_trj.gro' % name_ppm,
                                         skip=10, pbc_nojump=True, silent=True)
                        result = self.ppm_is_converged('%s_trj.gro' % name_ppm)
                        info_dict['failed'].append(result.get('failed'))
                        info_dict['continue'].append(result.get('continue'))
                        info_dict['continue_n'].append(result.get('continue_n'))
                        info_dict['reason'].append(result.get('reason'))
                        if more_info:
                            info_dict['more_info'].append(result.get('more_info'))
                        os.remove('%s_trj.gro' % name_ppm)
            ###

            inv_series = df['1/Viscosity']

            # select last 4/5 of data
            when = inv_series.index[len(inv_series) // 5]

            # use block average to estimate stderr, because 1/viscosity fluctuate heavily
            inv_blocks = average_of_blocks(inv_series.loc[when:])
            vis_blocks = [1000 / inv for inv in inv_blocks]  # convert Pa*s to cP
            a_list.append(ppm)
            vis_and_stderr = [np.mean(vis_blocks), np.std(vis_blocks, ddof=1) / math.sqrt(len(vis_blocks))]
            vis_list.append(vis_and_stderr[0])
            stderr_list.append(vis_and_stderr[1])
            if info_dict.get('failed')[-1]==False and info_dict.get('continue')[-1]==False and vis_and_stderr[1]/vis_and_stderr[0]>0.1:
                info_dict['continue'][-1] = True
                info_dict['continue_n'][-1] = int(1e7)
                info_dict['reason'][-1] = 'error bar too large for viscosity calculation'
        # if set(info_dict.get('failed'))=={False} and set(info_dict.get('continue'))=={False}:
        # coef_, score = polyfit(self.amplitudes_steps.keys(), vis_list, 1, weight=1 / np.sqrt(stderr_list))

        coef_, score = polyfit(a_list, vis_list, 1)
        c1, s1 = polyfit(a_list, (np.array(vis_list) + np.array(stderr_list)).tolist(), 1)
        c2, s1 = polyfit(a_list, (np.array(vis_list) - np.array(stderr_list)).tolist(), 1)
        ad_dict = {
            'viscosity': coef_[0],
            'vis-stderr': (c1[0] - c2[0]) / 2,
            'score': score,
            'vis_list': vis_list,
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

    @staticmethod
    def post_process(T_list, P_list, result_list, n_mol_list, **kwargs) -> (dict, str):
        t_set = set(T_list)
        p_set = set(P_list)
        def round3(x):
            return float('%.3e' % x)

        if len(p_set)==1:
            vis_score_list = [list(map(round3, [result['viscosity'], result['score']])) for result in result_list]
            vis_list = [i[0] for i in vis_score_list]
            t_p_vis_score_list = list(map(list, zip(T_list, P_list, vis_score_list)))
            t_p_vis_score_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            t_p_vis_list = list(map(list, zip(T_list, P_list, vis_list)))
            t_p_vis_list.sort(key=lambda x: (x[1], x[0]))  # sorted by P, then T
            p_coeff_score_list = []
            from ...analyzer.fitting import fit_VTF
            import numpy as np
            _t_list = [element[0] for element in t_p_vis_list]
            _vis_list = [element[2] for element in t_p_vis_list]
            _y_fit = np.log(_vis_list)
            _t_vis_coeff, _t_vis_score = fit_VTF(_t_list, _y_fit)
            # p_coeff_score_list.append([p, coeff, score])

            p = str(P_list[0])
            t_vis_VFT = {}
            t_min = _t_vis_coeff[1]
            t_max = 1000
            t_vis_VFT[p] = [list(map(round3, _t_vis_coeff)), round3(_t_vis_score), t_min, t_max]

            post_result = {
                't_p_vis_score_list': t_p_vis_score_list,
                'vis-t-VFT' : t_vis_VFT, # vis = np.exp(coeff[0]+coeff[2]/(t-coeff[1]))
                # 'p_coeff_score_list': p_coeff_score_list  # vis = np.exp(coeff[0]+coeff[2]/(t-coeff[1]))
            }
            return post_result, 'not important'

    @staticmethod
    def get_post_data(post_result, T, P, smiles_list, **kwargs) -> dict:
        import numpy as np
        converge_criterion = 0.95  # R value of fitting
        if len(post_result.get('vis-t-VFT'))==1:
            coef, score, tmin, tmax = post_result['vis-t-VFT'][str(P)]
            if not (score < converge_criterion or T < tmin - 10 or T > tmax + 10):
                vis = np.exp(coef[0]+coef[2]/(T-coef[1]))
            return {
                'viscosity': vis
            }
        else:
            return None



