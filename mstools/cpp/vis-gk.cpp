//#include "../../code_headfile/xy.h"
#include "/home/xiangyan/md1400/Github/cppheadfile/xy.h"
using namespace std;

int main(int argc, char** argv)
{
	FILE* f_pre= fopen(argv[1], "r");
	vector <double> t;
	vector <double> pxy;
	vector <double> pxz;
	vector <double> pyz;
	char p[1000];
	string split_symbol = " \n\t";
	while (fgets(p, 1000, f_pre)) {
		vector <string> sp = split(string(p), split_symbol);
		if (sp.size() <= 1) {

		}
		else if (sp[1] == "s0" or sp[1] == "s1" or sp[1] == "s2") {
			if (sp[3] != "\"Pres-XY\"" and sp[3] != "\"Pres-XZ\"" and sp[3] != "\"Pres-YZ\"") {
				error("input file error\n");
			}
		}
		else if (p[0] != '#' and p[0] != '@') {
			t.push_back(atof(sp[0].c_str()));
			pxy.push_back(atof(sp[1].c_str()));
			pxz.push_back(atof(sp[2].c_str()));
			pyz.push_back(atof(sp[3].c_str()));
		}
	}
	fclose(f_pre);
	double dt = t[1] - t[0];
	vector <double> t_list;
	vector <double> acf_list;
	for (unsigned i = 0; i < (t.size() / 2); ++i) {
		//printf("\r%i / %i", i, (t.size() / 2));
		double Dt = i * dt;
		double acf = 0.;
		for (unsigned j = 0; j < t.size() - i; ++j) {
			acf += pxy[j] * pxy[i + j];
			acf += pxz[j] * pxz[i + j];
			acf += pyz[j] * pyz[i + j];
		}
		acf /= (3. * (t.size() - i));
		t_list.push_back(Dt);
		acf_list.push_back(acf);
	}
	FILE* fout = fopen("acf.txt", "w");
	fprintf(fout, "#time(ps)\tACF(Pab)\n");
	for (unsigned i = 0; i < t_list.size(); ++i) {
		fprintf(fout, "%f\t%f\n", t_list[i], acf_list[i]);
	}
	fout = fopen("vis.txt", "w");
	fprintf(fout, "#time(ps)\tviscosity(mPa·s)\n");
	double V = atof(argv[2]);
	double T = atof(argv[3]);
	double convert = 6.022 * 0.001 * V / (8.314 * T);
	double vis = convert * acf_list[0] * dt / 2;
	for (unsigned i = 1; i < t_list.size(); ++i) {
		fprintf(fout, "%f\t%f\n", t_list[i] - 0.5 * dt, vis);
		vis += convert * acf_list[i] * dt;
	}
	fclose(fout);
}
