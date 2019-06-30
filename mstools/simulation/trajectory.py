mass_dict = {
    'H': 1.008,
    'B': 10.81,
    'C': 12.01,
    'N': 14.01,
    'O': 16.00,
    'F': 19.00,
    'Si':28.08,
    'P': 30.97,
    'S': 32.06,
    'Cl':35.45,
    'Br':79.90
}

class Frame():
    def __init__(self, atom_number = 0, boxX = 0., boxY = 0., boxZ = 0., t = 0., step = 0):
        self.atom_number = atom_number
        self.boxX = boxX
        self.boxY = boxY
        self.boxZ = boxZ
        self.t = t
        self.step = step
        self.atom_id = []
        self.mol_id = []
        self.atom_name = []
        self.mol_name = []
        self.x = []
        self.y = []
        self.z = []
        self.mass = []
        self.flag = []

    def check(self, mass_check=False):
        if self.atom_number < 1:
            raise Exception('Frame check error 1')
        if not (self.atom_number == len(self.atom_id) == len(self.atom_name) == len(self.mol_id) == len(self.mol_name) == len(self.x) == len(self.y) == len(self.z)):
            raise Exception('Frame check error 2')
        if self.boxX < 0. or self.boxY < 0. or self.boxZ < 0.:
            raise Exception('Frame check error 3')
        if mass_check:
            if not self.atom_number == len(self.mass):
                raise Exception('Frame check error 4')
    #
    def COMTransform(self):
        result = Frame(atom_number=self.atom_number, boxX=self.boxX, boxY=self.boxY, boxZ=self.boxZ, t=self.t, step=self.step)
        mid = 1
        mass = 0.
        x = y = z =0.
        for i in range(self.atom_number):
            if self.mol_id[i] == mid:
                x += self.x[i] * self.mass[i]
                y += self.y[i] * self.mass[i]
                z += self.z[i] * self.mass[i]
                mass += self.mass[i]
            elif self.mol_id[i] == mid+1:
                result.x.append(x / mass)
                result.y.append(y / mass)
                result.z.append(z / mass)
                result.mass.append(mass)
                result.atom_id.append(mid)
                result.atom_name.append(self.mol_name[i - 1])
                result.mol_id.append(mid)
                result.mol_name.append(self.mol_name[i - 1])
                mid += 1
                x = self.x[i] * self.mass[i]
                y = self.y[i] * self.mass[i]
                z = self.z[i] * self.mass[i]
                mass = self.mass[i]
            if i == self.atom_number - 1:
                result.x.append(x / mass)
                result.y.append(y / mass)
                result.z.append(z / mass)
                result.mass.append(mass)
                result.atom_id.append(mid)
                result.atom_name.append(self.mol_name[i])
                result.mol_id.append(mid)
                result.mol_name.append(self.mol_name[i])
                break
        result.atom_number = len(result.x)
        self.atom_number = result.atom_number
        self.x = result.x
        self.y = result.y
        self.z = result.z
        self.mass = result.mass
        self.atom_id = result.atom_id
        self.atom_name = result.atom_name
        self.mol_id = result.mol_id
        self.mol_name = result.mol_name

    # present the frame in center of atoms
    def CenterTransform(self, mass_flag=False):
        if mass_flag:
            if None in self.mass:
                raise Exception('The atom is not in mass_dict')
        result = Frame(atom_number=self.atom_number, boxX=self.boxX, boxY=self.boxY, boxZ=self.boxZ, t=self.t, step=self.step)
        mid = 1
        atom_n_in_mol = 0
        if mass_flag:
            mass = 0.
        x = y = z =0.
        for i in range(self.atom_number):
            if self.mol_id[i] == mid:
                atom_n_in_mol += 1
                x += self.x[i]
                y += self.y[i]
                z += self.z[i]
                if mass_flag:
                    mass += self.mass[i]
            elif self.mol_id[i] == mid+1:
                result.x.append(x / atom_n_in_mol)
                result.y.append(y / atom_n_in_mol)
                result.z.append(z / atom_n_in_mol)
                if mass_flag:
                    result.mass.append(mass)
                result.atom_id.append(mid)
                result.atom_name.append(self.mol_name[i - 1])
                result.mol_id.append(mid)
                result.mol_name.append(self.mol_name[i - 1])
                mid += 1
                x = self.x[i]
                y = self.y[i]
                z = self.z[i]
                if mass_flag:
                    mass = self.mass[i]
                atom_n_in_mol = 1
            if i == self.atom_number - 1:
                result.x.append(x / atom_n_in_mol)
                result.y.append(y / atom_n_in_mol)
                result.z.append(z / atom_n_in_mol)
                if mass_flag:
                    result.mass.append(mass)
                result.atom_id.append(mid)
                result.atom_name.append(self.mol_name[i])
                result.mol_id.append(mid)
                result.mol_name.append(self.mol_name[i])
                break
        result.atom_number = len(result.x)
        self.atom_number = result.atom_number
        self.x = result.x
        self.y = result.y
        self.z = result.z
        self.mass = result.mass
        self.atom_id = result.atom_id
        self.atom_name = result.atom_name
        self.mol_id = result.mol_id
        self.mol_name = result.mol_name

    # transform the coordinates in reduced units in box size
    def ReducedUnitTransform(self):
        for i in range(self.atom_number):
            self.x[i] /= self.boxX
            self.y[i] /= self.boxY
            self.z[i] /= self.boxZ

class Trajectory():
    def __init__(self, traj_name, readmass=False, COM=False, center=False, silent=True):
        self.traj_name = traj_name
        self.readmass = readmass
        self.COM = COM
        self.center = center
        self.silent = silent
        self.traj_info = []
        self.use_check()
        self.ReadTraj()

    def use_check(self):
        if self.COM == True and self.readmass == False:
            raise Exception('you need to set readmass=True, for COM=True\n')


    def ReadTraj(self):
        if self.traj_name.split('.')[-1]=='gro':
            self.ReadGRO()

    def ReadGRO(self):
        lines = open(self.traj_name).readlines()
        atom_number = int(lines[1])
        frame_number = int(len(lines) / (atom_number+3))
        for i in range(frame_number):
            if not self.silent:
                print('reading %i-th frame' % (i+1))
            temp_frame = Frame()
            line_info = lines[i * (atom_number+3)].split()
            for j, info in enumerate(line_info):
                if info =='t=':
                    temp_frame.t = float(line_info[j+1])
                elif info =='step=':
                    temp_frame.step = int(line_info[j+1])

            temp_frame.atom_number = int( lines[i * (atom_number+3) + 1] )
            temp_frame.boxX = float(lines[(i + 1) * (atom_number + 3) - 1].split()[0])
            temp_frame.boxY = float(lines[(i + 1) * (atom_number + 3) - 1].split()[1])
            temp_frame.boxZ = float(lines[(i + 1) * (atom_number + 3) - 1].split()[2])
            for j in range(atom_number):
                k = i * (atom_number+3) + 2 +j
                temp = list(lines[k])
                temp.insert(15, ' ')
                temp.insert(5, ' ')
                temp = "".join(temp).split()
                temp_frame.mol_id.append(int(temp[0]))
                temp_frame.mol_name.append(temp[1])
                temp_frame.atom_name.append(temp[2])
                temp_frame.atom_id.append(int(temp[3]))
                temp_frame.x.append(float(temp[4]))
                temp_frame.y.append(float(temp[5]))
                temp_frame.z.append(float(temp[6]))
                if self.readmass:
                    mass = mass_dict.get(temp[2])
                    if mass == None:
                        raise Exception('expect mass info for atom: %s\n' % (temp[2]))
                    else:
                        temp_frame.mass.append(mass)
            if self.COM:
                temp_frame.COMTransform()
            elif self.center:
                temp_frame.CenterTransform(mass_flag=self.readmass)
            self.traj_info.append(temp_frame)
            self.traj_info[-1].check(mass_check=self.readmass)

    def ReducedUnitTransform(self):
        for frame in self.traj_info:
            frame.ReducedUnitTransform()




