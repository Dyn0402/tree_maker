/*
 * BES_Pars.h
 *
 *  Created on: Jul 16, 2020
 *      Author: Dylan Neff
 */

#ifndef BESPARS_H_
#define BESPARS_H_

#include <map>
#include <vector>

using namespace std;


class BESPars {
public:
	// Structors
	BESPars();
	BESPars(int energy);
	BESPars(int energy, int bes_phase);
	~BESPars();

	// Getters

	// Setters
	void set_energy(int energy);
	void set_bes_phase(int bes_phase);
	void set_energy_bes(int energy, int bes_phase);

	// Attributes
	int branch_buffer = 500000;
	int branch_split = 1;

	float vertex_min = 1e-5;

	int min_14GeV_run = 15049041;  // From Toshihiro's code

	float vz_max;
	float vr_max;
	float vy_offset;
	float vpd_vz_max_diff;
	vector<int> triggers;
	vector<int> bad_runs;


private:
	// Attributes
	int energy;
	int bes_phase;

	// Doers
	void set_pars();

	// Parameter Values  All QA taken from Analysis note for net-proton C4 fluctuations 2019

	map<int, map<int, float>> vz_max_map {  // [energy][bes_phase]
		{7, {{1, 50.0}, {2, 50.0}}},
		{11, {{1, 30.0}, {2, 30.0}}},
		{15, {{1, 30.0}, {2, 30.0}}},
		{19, {{1, 30.0}, {2, 30.0}}},
		{27, {{1, 30.0}, {2, 30.0}}},
		{39, {{1, 30.0}, {2, 30.0}}},
		{54, {{1, 30.0}, {2, 30.0}}},
		{62, {{1, 30.0}, {2, 30.0}}},
		{200, {{1, 30.0}, {2, 30.0}}}
	};

	map<int, map<int, float>> vr_max_map {
		{7, {{1, 2.0}, {2, 2.0}}},
		{11, {{1, 2.0}, {2, 2.0}}},
		{15, {{1, 1.0}, {2, 1.0}}},
		{19, {{1, 2.0}, {2, 2.0}}},
		{27, {{1, 2.0}, {2, 2.0}}},
		{39, {{1, 2.0}, {2, 2.0}}},
		{54, {{1, 2.0}, {2, 2.0}}},
		{62, {{1, 2.0}, {2, 2.0}}},
		{200, {{1, 2.0}, {2, 2.0}}}
	};

	map<int, map<int, float>> vy_offset_map {
		{7, {{1, 0.0}, {2, 0.0}}},
		{11, {{1, 0.0}, {2, 0.0}}},
		{15, {{1, 0.89}, {2, 0.89}}},
		{19, {{1, 0.0}, {2, 0.0}}},
		{27, {{1, 0.0}, {2, 0.0}}},
		{39, {{1, 0.0}, {2, 0.0}}},
		{54, {{1, 0.0}, {2, 0.0}}},
		{62, {{1, 0.0}, {2, 0.0}}},
		{200, {{1, 0.0}, {2, 0.0}}}
	};

	map<int, map<int, float>> vpd_vz_max_diff_map {
//		{7, {{1, -1.0}, {2, -1.0}}},  // - sign in initialization gives compile bad char in name error. Think it's a ROOT issue, very worrisome.
//		{11, {{1, -1.0}, {2, -1.0}}},  // -1 ignore code
//		{15, {{1, -1.0}, {2, -1.0}}},
//		{19, {{1, -1.0}, {2, -1.0}}},
//		{27, {{1, -1.0}, {2, -1.0}}},
		{39, {{1, 3.0}, {2, 3.0}}},
		{54, {{1, 3.0}, {2, 3.0}}},
		{62, {{1, 3.0}, {2, 3.0}}},
		{200, {{1, 3.0}, {2, 3.0}}}
	};

	map<int, map<int, vector<int>>> triggers_map {
		{7, {
				{1, {290001, 290004}},  // Conflict with Roli Original
				{2, {}}}},
		{11, {
				{1, {310004, 310014}},  // Conflict with Xiaofeng Analysis Note
				{2, {}}}},
		{15, {
				{1, {440005, 440006, 440015, 440016}},  // Conflict with Roli Original
				{2, {}}}},
		{19, {
				{1, {340001, 340011, 340021}},
				{2, {}}}},
		{27, {
				{1, {360001}},
				{2, {}}}},
		{39, {
				{1, {280001}},
				{2, {}}}},
		{54, {
				{1, {580001, 580021}},  // From Xiaofeng Analysis Note, Roli didn't have 54
				{2, {}}}},
		{62, {
				{1, {270001, 270011, 270021}},
				{2, {}}}},
		{200, {
				{1, {350001, 350011, 350003, 350013, 350023, 350033, 350043}},  // Conflict with Xiaofeng Analysis Note
				{2, {}}}}
	};

	// Compiled from Xiaofeng's student Yu, still some runs in StRefMultCorr that aren't here so not final.
	map<int, map<int, vector<int>>> bad_runs_map {
		{7, {
				{1, {11114074, 11114084, 11114085, 11114086, 11114088, 11114089, 11114091, 11114094, 11114095, 11114096, 11114098, 11114100, 11114109, 11114111, 11115005, 11115007, 11115010, 11115011, 11115013, 11115014, 11115019, 11115020, 11115025, 11115026, 11115027, 11115028, 11115030, 11115031, 11115032, 11115033, 11115036, 11115042, 11115044, 11115050, 11115051, 11115062, 11115064, 11115068, 11115069, 11115072, 11115078, 11115079, 11115080, 11115084, 11115086, 11115088, 11115093, 11115094, 11116001, 11116002, 11116005, 11116006, 11116010, 11116014, 11116020, 11116023, 11116027, 11116028, 11116060, 11116061, 11116062, 11116063, 11116064, 11116068, 11116070, 11116071, 11116072, 11116073, 11116075, 11117002, 11117006, 11117020, 11117026, 11117029, 11117030, 11117031, 11117032, 11117033, 11117034, 11117036, 11117038, 11117039, 11117041, 11117043, 11117044, 11117045, 11117046, 11117052, 11117054, 11117055, 11117063, 11117064, 11117071, 11117072, 11117074, 11117075, 11117076, 11117085, 11117088, 11117089, 11117090, 11117093, 11117094, 11117095, 11117098, 11117100, 11117102, 11117103, 11117104, 11117107, 11118007, 11118008, 11118016, 11118023, 11118024, 11118025, 11118026, 11118039, 11118043, 11118044, 11119001, 11119003, 11119006, 11119007, 11119009, 11119012, 11119013, 11119015, 11119016, 11119017, 11119020, 11119022, 11119024, 11119025, 11119026, 11119029, 11119030, 11119054, 11119056, 11119057, 11119060, 11119061, 11119062, 11119064, 11119066, 11119067, 11119068, 11119069, 11119070, 11119071, 11119074, 11119075, 11119077, 11119079, 11119081, 11119082, 11119087, 11119090, 11119091, 11119094, 11119100, 11119101, 11120003, 11120006, 11120008, 11120011, 11120014, 11120019, 11120022, 11120023, 11120028, 11120030, 11120034, 11120037, 11120039, 11120040, 11120042, 11120043, 11120045, 11120049, 11120050, 11120052, 11120057, 11120060, 11120062, 11120063, 11120069, 11120070, 11120071, 11120074, 11120077, 11120078, 11120083, 11120084, 11120089, 11120090, 11120092, 11121002, 11121006, 11121013, 11121014, 11121015, 11121019, 11121029, 11121030, 11121031, 11121034, 11121035, 11121043, 11121044, 11121049, 11121054, 11121058, 11121065, 11121066, 11121067, 11121070, 11121075, 11121082, 11122001, 11122007, 11122008, 11122010, 11122017, 11122024, 11122037, 11122038, 11122047, 11122048, 11122049, 11122050, 11122051, 11122052, 11122053, 11122058, 11122061, 11122062, 11122067, 11122068, 11122069, 11122073, 11122074, 11122078, 11122079, 11122082, 11122085, 11122097, 11123003, 11123004, 11123015, 11123026, 11123028, 11123040, 11123044, 11123051, 11123053, 11123054, 11123055, 11123056, 11123057, 11123058, 11123059, 11123060, 11123067, 11123075, 11123076, 11123077, 11123079, 11123081, 11123082, 11123084, 11123086, 11123088, 11123089, 11123090, 11123091, 11123093, 11123094, 11123095, 11123098, 11123100, 11123101, 11123102, 11123104, 11123106, 11124001, 11124002, 11124005, 11124006, 11124007, 11124008, 11124010, 11124011, 11124014, 11124015, 11124016, 11124018, 11124019, 11124037, 11124041, 11124044, 11124046, 11124047, 11124048, 11124049, 11124050, 11124051, 11124052, 11124053, 11124054, 11124055, 11124058, 11124060, 11124061, 11124062, 11124063, 11124064, 11124065, 11124066, 11124067, 11124069, 11125001, 11125002, 11125003, 11125004, 11125005, 11125006, 11125007, 11125008, 11125012, 11125013, 11125014, 11125015, 11125016, 11125017, 11125019, 11125020, 11125021, 11125022, 11125023, 11125024, 11125073, 11125081, 11125089, 11125090, 11125093, 11125096, 11125097, 11126005, 11126006, 11126007, 11126016, 11126017, 11126018, 11126022, 11126023, 11126024, 11127001, 11127002, 11127020, 11127039, 11127043, 11128005, 11128012, 11128018, 11128045, 11128050, 11128056, 11128072, 11129011, 11129018, 11129022, 11129025, 11129028, 11129051, 11130027, 11130034, 11130057, 11130080, 11131038, 11131062, 11132013, 11132070, 11133006, 11133016, 11133019, 11134032, 11134053, 11134060, 11134067, 11134076, 11135030, 11135059, 11135068, 11136003, 11136005, 11136006, 11136007, 11136008, 11136012, 11136013, 11136014, 11136015, 11136019, 11136030, 11136061, 11136067, 11136069, 11136073, 11136076, 11136101, 11136130, 11136160, 11136163, 11137019, 11137021, 11137044, 11137049, 11137116, 11138012, 11138016, 11138022, 11138027, 11138049, 11138069, 11138086, 11138124, 11139014, 11140076, 11140086, 11141063, 11142091, 11142117, 11143016, 11143026, 11143028, 11144001, 11144009, 11144031, 11144033, 11144040, 11144043, 11144052, 11145008, 11145028, 11145035, 11146061, 11146076, 11146079, 11146095, 11147004, 11147005, 11147006, 11147014, 11147017, 11147021, 11147023, 11114074, 11115017, 11115066, 11116011, 11117003, 11117058, 11118001, 11119020, 11119083, 11120021, 11120068, 11121013, 11121045, 11122012, 11122065, 11122111, 11123047, 11123084, 11124009, 11125008, 11129025, 11136076, 11114076, 11115018, 11115068, 11116012, 11117004, 11117060, 11118002, 11119021, 11119085, 11120022, 11120069, 11121014, 11121046, 11122015, 11122066, 11122112, 11123048, 11123085, 11124010, 11125012, 11129027, 11136130, 11114084, 11115019, 11115069, 11116013, 11117005, 11117061, 11118006, 11119022, 11119087, 11120023, 11120070, 11121015, 11121048, 11122016, 11122067, 11122113, 11123050, 11123086, 11124011, 11125013, 11129028, 11136163, 11114085, 11115020, 11115070, 11116014, 11117006, 11117062, 11118007, 11119024, 11119088, 11120027, 11120071, 11121016, 11121049, 11122017, 11122068, 11123001, 11123051, 11123087, 11124012, 11125014, 11129031, 11137019, 11114086, 11115021, 11115071, 11116015, 11117007, 11117063, 11118008, 11119025, 11119090, 11120028, 11120072, 11121019, 11121050, 11122018, 11122069, 11123003, 11123053, 11123088, 11124013, 11125015, 11129033, 11138049, 11114087, 11115022, 11115072, 11116016, 11117009, 11117064, 11118015, 11119026, 11119091, 11120030, 11120073, 11121021, 11121054, 11122019, 11122070, 11123004, 11123054, 11123089, 11124014, 11125016, 11129036, 11138086, 11114088, 11115023, 11115073, 11116017, 11117020, 11117067, 11118016, 11119029, 11119093, 11120031, 11120074, 11121022, 11121057, 11122020, 11122071, 11123005, 11123055, 11123090, 11124015, 11125017, 11129038, 11138124, 11114089, 11115024, 11115075, 11116019, 11117021, 11117068, 11118023, 11119030, 11119094, 11120034, 11120077, 11121023, 11121058, 11122021, 11122073, 11123006, 11123056, 11123091, 11124016, 11125019, 11129044, 11140012, 11114091, 11115025, 11115076, 11116020, 11117026, 11117069, 11118024, 11119033, 11119098, 11120035, 11120078, 11121024, 11121059, 11122023, 11122075, 11123007, 11123057, 11123092, 11124017, 11125020, 11129045, 11140076, 11114094, 11115026, 11115078, 11116021, 11117027, 11117070, 11118025, 11119054, 11119099, 11120036, 11120079, 11121025, 11121060, 11122024, 11122077, 11123010, 11123058, 11123093, 11124018, 11125021, 11129046, 11140088, 11114095, 11115027, 11115079, 11116023, 11117028, 11117071, 11118026, 11119055, 11119100, 11120037, 11120081, 11121026, 11121061, 11122025, 11122078, 11123011, 11123059, 11123094, 11124019, 11125022, 11129047, 11140094, 11114096, 11115028, 11115080, 11116027, 11117029, 11117072, 11118039, 11119056, 11119101, 11120038, 11120083, 11121027, 11121062, 11122037, 11122079, 11123012, 11123060, 11123095, 11124035, 11125023, 11129049, 11141053, 11114097, 11115030, 11115082, 11116028, 11117030, 11117073, 11118040, 11119057, 11119102, 11120039, 11120084, 11121028, 11121063, 11122038, 11122080, 11123013, 11123062, 11123097, 11124037, 11125024, 11129050, 11141063, 11114098, 11115031, 11115084, 11116031, 11117031, 11117074, 11118043, 11119060, 11120002, 11120040, 11120085, 11121029, 11121064, 11122039, 11122081, 11123014, 11123063, 11123098, 11124041, 11125070, 11129051, 11142117, 11114100, 11115032, 11115085, 11116054, 11117032, 11117075, 11118044, 11119061, 11120003, 11120041, 11120086, 11121030, 11121066, 11122040, 11122082, 11123015, 11123064, 11123099, 11124044, 11125071, 11129061, 11143016, 11114109, 11115033, 11115086, 11116055, 11117033, 11117076, 11119001, 11119062, 11120004, 11120042, 11120087, 11121031, 11121067, 11122047, 11122083, 11123016, 11123065, 11123100, 11124046, 11125072, 11129063, 11143026, 11114111, 11115034, 11115087, 11116056, 11117034, 11117079, 11119003, 11119063, 11120006, 11120043, 11120088, 11121034, 11121068, 11122048, 11122084, 11123017, 11123066, 11123101, 11124047, 11125073, 11129067, 11143028, 11114112, 11115036, 11115088, 11116057, 11117035, 11117085, 11119004, 11119064, 11120007, 11120045, 11120089, 11121042, 11121070, 11122049, 11122085, 11123020, 11123067, 11123102, 11124048, 11125075, 11129075, 11144001, 11115001, 11115039, 11115089, 11116060, 11117036, 11117088, 11119006, 11119066, 11120008, 11120046, 11120090, 11121043, 11121071, 11122050, 11122086, 11123021, 11123068, 11123103, 11124049, 11125076, 11129081, 11144031, 11115004, 11115041, 11115093, 11116061, 11117037, 11117089, 11119007, 11119067, 11120009, 11120047, 11120091, 11121044, 11121072, 11122051, 11122087, 11123022, 11123069, 11123104, 11136007, 11125078, 11130002, 11144033, 11115005, 11115042, 11115094, 11116062, 11117038, 11117090, 11119009, 11119068, 11120010, 11120050, 11120092, 11121074, 11121088, 11122052, 11122088, 11123023, 11123071, 11123105, 11136008, 11125081, 11130003, 11144043, 11115006, 11115043, 11116001, 11116063, 11117039, 11117093, 11119010, 11119070, 11120011, 11120051, 11121001, 11121075, 11121089, 11122053, 11122089, 11123024, 11123072, 11123106, 11136012, 11125085, 11130004, 11145008, 11115007, 11115044, 11116002, 11116064, 11117041, 11117094, 11119011, 11119071, 11120012, 11120052, 11121002, 11121076, 11121090, 11122054, 11122090, 11123025, 11123075, 11123107, 11136013, 11125089, 11130005, 11145028, 11115008, 11115050, 11116003, 11116068, 11117043, 11117095, 11119012, 11119073, 11120013, 11120057, 11121003, 11121077, 11122001, 11122055, 11122092, 11123026, 11123076, 11124001, 11136014, 11125090, 11130013, 11145035, 11115010, 11115051, 11116004, 11116070, 11117044, 11117098, 11119013, 11119074, 11120014, 11120058, 11121004, 11121078, 11122005, 11122057, 11122094, 11123027, 11123077, 11124002, 11136015, 11125093, 11130014, 11145042, 11115011, 11115052, 11116005, 11116071, 11117045, 11117100, 11119014, 11119075, 11120015, 11120059, 11121005, 11121080, 11122006, 11122058, 11122095, 11123028, 11123078, 11124003, 11136020, 11125094, 11130015, 11146059, 11115012, 11115053, 11116006, 11116072, 11117046, 11117102, 11119015, 11119076, 11120016, 11120060, 11121006, 11121081, 11122007, 11122059, 11122096, 11123029, 11123079, 11124004, 11136030, 11125096, 11130017, 11147004, 11115013, 11115061, 11116007, 11116073, 11117052, 11117103, 11119016, 11119077, 11120017, 11120062, 11121007, 11121082, 11122008, 11122060, 11122097, 11123039, 11123080, 11124005, 11136069, 11125097, 11130022, 11147005, 11115014, 11115062, 11116008, 11116074, 11117054, 11117104, 11119017, 11119079, 11120018, 11120063, 11121008, 11121083, 11122009, 11122061, 11122098, 11123040, 11123081, 11124006, 11126012, 11130028, 11125099, 11130023, 11115015, 11115064, 11116009, 11116075, 11117055, 11117106, 11119018, 11119081, 11120019, 11120066, 11121011, 11121084, 11122010, 11122062, 11122099, 11123044, 11123082, 11124007, 11126016, 11130048, 11126005, 11130024, 11115016, 11115065, 11116010, 11116076, 11117057, 11117107, 11119019, 11119082, 11120020, 11120067, 11121012, 11121086, 11122011, 11122064, 11122101, 11123046, 11123083, 11124008, 11126006, 11130025, 11147017, 11147006, 11124050, 11126017, 11130050, 11124067, 11128014, 11124051, 11126018, 11130057, 11124069, 11128045, 11124052, 11126022, 11131004, 11125001, 11128050, 11124053, 11126023, 11131038, 11125002, 11128056, 11147014, 11136005, 11126007, 11130027, 11147021, 11124055, 11127001, 11132070, 11125004, 11129010, 11124060, 11127028, 11133016, 11125006, 11129018, 11124061, 11127033, 11134053, 11125007, 11129022, 11124062, 11127039, 11134061, 11124065, 11128005, 11136003, 11128012, 11124058, 11127002, 11133006, 11125005, 11129011, 11124054, 11126044, 11131062, 11125003, 11129001, 11124063, 11127043, 11134067, 11124064, 11127054, 11135068, 11124066}},
				{2, {}}}},
		{11, {
				{1, {11148039, 11148045, 11148055, 11148069, 11149001, 11149008, 11149010, 11149011, 11149015, 11149017, 11149018, 11149047, 11150016, 11150017, 11150025, 11150028, 11150029, 11151036, 11151040, 11151050, 11151057, 11151084, 11151086, 11152016, 11152036, 11152078, 11153032, 11153042, 11153045, 11154026, 11154034, 11155001, 11155009, 11156003, 11156008, 11156009, 11156036, 11156043, 11156044, 11156045, 11157012, 11158006, 11158022, 11158024}},
				{2, {}}}},
		{15, {
				{1, {15053027, 15056125, 15062006, 15054053, 15057014, 15053028, 15057001, 15063017, 15054054, 15057015, 15053029, 15057002, 15065014, 15055131, 15057016, 15053034, 15057003, 15066008, 15055133, 15057018, 15053048, 15057004, 15066013, 15055134, 15057019, 15053049, 15057005, 15066017, 15056113, 15057055, 15053050, 15057006, 15066018, 15056114, 15060061, 15053052, 15057007, 15066071, 15056116, 15060062, 15053053, 15057008, 15066072, 15056117, 15061001, 15053054, 15057010, 15066074, 15056119, 15061002, 15053056, 15057011, 15066082, 15056123, 15061009, 15053057, 15057012, 15067027, 15056124, 15062002, 15053058, 15057013, 15070009}},
				{2, {}}}},
		{19, {
				{1, {12113084, 12113091, 12114007, 12114063, 12114092, 12114101, 12114116, 12115009, 12115014, 12115015, 12115016, 12115018, 12115019, 12115020, 12115022, 12115023, 12115073, 12115093, 12116012, 12116054, 12117010, 12117016, 12117047, 12118036, 12119023, 12119032, 12119039, 12119040, 12119042, 12120017, 12120018, 12120026, 12121017, 12121022, 12121034, 12122019}},
				{2, {}}}},
		{27, {
				{1, {12172049, 12172050, 12172051, 12172055, 12173009, 12173018, 12173026, 12173030, 12173031, 12173032, 12173033, 12173034, 12173047, 12173054, 12173055, 12173056, 12174076, 12174085, 12174086, 12175007, 12175030, 12175062, 12175087, 12175113, 12175114, 12175115, 12176001, 12176044, 12176046, 12176047, 12176054, 12176067, 12176104, 12177015, 12177061, 12177092, 12177097, 12177099, 12177101, 12177106, 12177107, 12177108, 12177126, 12178003, 12178004, 12178005, 12178006, 12178051, 12178093, 12178099, 12178120, 12179085, 12179094}},
				{2, {}}}},
		{39, {
				{1, {11099124, 11100002, 11100045, 11101046, 11102012, 11102051, 11102052, 11102053, 11102054, 11102055, 11102058, 11102098, 11103007, 11103008, 11103009, 11103035, 11103047, 11103056, 11103058, 11103065, 11103092, 11103093, 11105052, 11105053, 11105054, 11105055, 11106026, 11106027, 11106028, 11106039, 11106040, 11106041, 11107007, 11107042, 11107057, 11107061, 11107065, 11107074, 11108004, 11108040, 11108053, 11108065, 11108101, 11109013, 11109077, 11109088, 11109090, 11109102, 11109105, 11109127, 11110005, 11110013, 11110034, 11110042, 11110073, 11110076, 11111084, 11111085}},
				{2, {}}}},
		{54, {
				{1, {}},
				{2, {}}}},
		{62, {
				{1, {11080054, 11080055, 11080057, 11080060, 11080061, 11080062, 11080064, 11080070, 11080072, 11081001, 11081002, 11081003, 11081023, 11081025, 11081032, 11081033, 11081034, 11081052, 11082012, 11082013, 11082046, 11082056, 11082057, 11084009, 11084011, 11084012, 11084013, 11084020, 11084021, 11084035, 11084044, 11084064, 11085013, 11085015, 11085025, 11085030, 11085046, 11085055, 11085056, 11085057, 11086005, 11086007, 11087001, 11087002, 11087003, 11087004, 11087057, 11087058, 11087059, 11088013, 11089026, 11089028, 11089029, 11089048, 11089054, 11089055, 11089068, 11089072, 11091007, 11091015, 11091021, 11091078, 11092010, 11092011, 11092012, 11092031, 11092032, 11092033, 11092034, 11092067, 11092096, 11093001, 11094004, 11094016, 11094017, 11094018, 11094019, 11094020, 11094021, 11094022, 11094023, 11094024, 11094027, 11094028, 11094042, 11094044, 11094045, 11094046, 11094047, 11094048, 11094050, 11094051, 11094052, 11094053, 11094054, 11094055, 11094074, 11094075, 11094077, 11095001, 11095002, 11095003, 11095004, 11095005, 11095006, 11095009, 11095010, 11095011, 11095012, 11095013, 11095014, 11095015, 11095022, 11095040, 11095048, 11095049, 11095050, 11095051, 11095061, 11095062, 11095063, 11095064, 11095065, 11095082, 11095087, 11096024, 11096039, 11096043, 11096044, 11097093}},
				{2, {}}}},
		{200, {
				{1, {11002145, 11003001, 11003002, 11003004, 11003005, 11003006, 11003010, 11003011, 11003012, 11003013, 11003014, 11003015, 11003016, 11003017, 11003029, 11003030, 11003031, 11003062, 11003063, 11003064, 11003065, 11003066, 11003067, 11003068, 11003101, 11003102, 11004007, 11004008, 11004009, 11004010, 11004011, 11004012, 11004013, 11004014, 11004015, 11004016, 11004018, 11004020, 11004021, 11004023, 11004024, 11004025, 11004026, 11004028, 11004029, 11004030, 11004032, 11004033, 11004034, 11004035, 11004037, 11004038, 11004040, 11004041, 11005030, 11005042, 11006004, 11006005, 11006008, 11007006, 11007007, 11007008, 11007015, 11007019, 11009014, 11009032, 11010031, 11011019, 11011053, 11014041, 11014066, 11015069, 11015071, 11016024, 11016056, 11017006, 11018003, 11018007, 11018008, 11018036, 11019001, 11019077, 11019080, 11019081, 11021027, 11021028, 11021031, 11023021, 11023044, 11023046, 11023048, 11025032, 11025033, 11025034, 11025036, 11025037, 11025038, 11025039, 11025040, 11025051, 11025054, 11025067, 11025069, 11026005, 11026008, 11026010, 11026021, 11026022, 11026023, 11026025, 11026035, 11026041, 11026042, 11026067, 11026068, 11026114, 11028004, 11028005, 11028006, 11028007, 11028008, 11028009, 11028010, 11028011, 11028012, 11028013, 11028018, 11028019, 11028020, 11028021, 11028022, 11028023, 11028024, 11028025, 11028026, 11028027, 11030041, 11030042, 11030080, 11031022, 11031046, 11031047, 11031052, 11031061, 11031064, 11033039, 11035008, 11035009, 11035072, 11036026, 11037035, 11037037, 11037058, 11037060, 11037066, 11037067, 11038047, 11038048, 11038049, 11038050, 11039047, 11039067, 11040078, 11040083, 11041022, 11041023, 11041040, 11041041, 11042001, 11042002, 11042003, 11042004, 11042005, 11042006, 11042007, 11042008, 11042011, 11042012, 11042018, 11042019, 11042020, 11042021, 11042022, 11042023, 11042024, 11042025, 11042026, 11042027, 11042042, 11042043, 11042044, 11042045, 11042046, 11042047, 11042048, 11042049, 11044029, 11047059, 11047065, 11047066, 11047067, 11048037, 11049001, 11049002, 11049005, 11049023, 11051034, 11051037, 11051038, 11051049, 11051051, 11051055, 11051063, 11051064, 11051068, 11052011, 11053045, 11053057, 11054021, 11054022, 11054024, 11054058, 11054059, 11054062, 11054065, 11054066, 11055099, 11055102, 11056017, 11056029, 11056035, 11056040, 11057012, 11057019, 11057027, 11057030, 11057035, 11057036, 11057038, 11057041, 11057046, 11057048, 11057049, 11058002, 11058005, 11058015, 11058040, 11058050, 11058057, 11058058, 11058070, 11058073, 11058083, 11059006, 11059012, 11059013, 11059032, 11059043, 11059047, 11059055, 11059060, 11059061, 11059069, 11059075, 11059076, 11059077, 11060007, 11060008, 11060015, 11060034, 11060046, 11060047, 11060048, 11060049, 11060052, 11060053, 11060059, 11060069, 11060074, 11060075, 11060076, 11061008, 11061009, 11061014, 11061017, 11061021, 11061034, 11061037, 11061038, 11061039, 11061042, 11061043, 11061067, 11061080, 11061095, 11062011, 11062017, 11063006, 11063007, 11063008, 11063011, 11063013, 11063014, 11063015, 11063016, 11063017, 11063036, 11063060, 11063064, 11063070, 11063076, 11063077, 11063083, 11063086, 11064003, 11064007, 11064023, 11065038, 11066024, 11066045, 11067024, 11071056, 11072030, 11072031, 11072032, 11072044, 11072045, 11073001, 11073002, 11073003, 11073004, 11073049, 11075002, 11075004, 11075035, 11075039, 11075040, 11075045, 11075048, 11075053, 11075055, 11075063, 11077018}},
				{2, {}}}}
	};
};






#endif /* BESPARS_H_ */
