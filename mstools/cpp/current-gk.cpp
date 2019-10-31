//#include "../../code_headfile/xy.h"
#include "/home/xiangyan/md1400/Github/cppheadfile/xy.h"
using namespace std;

int main(int argc, char** argv)
{
	FILE* f_current= fopen(argv[1], "r");
	double V = atof(argv[2]);
	double T = atof(argv[3]);
	double weight = atof(argv[4]);
	char weight_char[10];
	sprintf(weight_char, "%.2f", weight);
	vector <double> t;
	vector <Vector> J;
	char p[1000];
	string split_symbol = " \n\t";
	while (fgets(p, 1000, f_current)) {
		vector <string> sp = split(string(p), split_symbol);
		if (sp.size() <= 1) {

		}
		else if (p[0] != '#' and p[0] != '@') {
			t.push_back(atof(sp[0].c_str()));
			Vector j(atof(sp[1].c_str()), atof(sp[2].c_str()), atof(sp[3].c_str()));
			J.push_back(j)
		}
	}
	fclose(f_current);
	double dt = t[1] - t[0];
	vector <double> t_list;
	vector <double> acf_list;
	for (unsigned i = 0; i < (t.size() / 4); ++i) {
		//printf("\r%i / %i", i, (t.size() / 2));
		double Dt = i * dt;
		double acf = 0.;
		for (unsigned j = 0; j < t.size() - i; ++j) {
			acf += dotProduct(J[j], J[i + j]);
		}
		acf /= (t.size() - i);
		t_list.push_back(Dt);
		acf_list.push_back(acf);
	}

	string fn = "J_acf.txt";
	FILE* fout = fopen(fn.c_str(), "w");
	fprintf(fout, "#time(ps)\tACF(J)\n");
	for (unsigned i = 0; i < t_list.size(); ++i) {
		fprintf(fout, "%f\t%f\n", t_list[i], acf_list[i]);
	}
	fn = "econ-" + string(weight_char) + ".txt";
	fout = fopen(fn.c_str(), "w");
	fprintf(fout, "#time(ps)\telectrical_conductivity(S/m)\n");

	double convert = 1.6 ** 2 * 6.022 * 10 ** 6 / (3 * 8.314 * T * V);
	double econ = convert * acf_list[0] * dt / 2;
	for (unsigned i = 1; i < t_list.size(); ++i) {
		fprintf(fout, "%f\t%f\n", t_list[i] - 0.5 * dt, econ);
		if(t_list[i] <= 1){
		    econ += convert * acf_list[i] * dt;
		}
		else{
		    econ += convert * acf_list[i] * dt * pow(t_list[i], -weight);
		}
	}
	fclose(fout);
}
