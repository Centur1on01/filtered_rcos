#define _USE_MATH_DEFINES

#include <iostream>
#include <complex>
#include <vector>
#include <deque>
#include <ctime>
#include <cmath>
#include <random>
#include <string>
#include <fstream>
#include <cstdint>

using namespace std;



class signalSource {
	// наследуемые классы: источник бита и маппер
	// источник бита берет из файла или m-последовательности
	// мапппер позволяет выбрать сигнальное созвездие

	// получение битов из файла -> выбор созвездия -> перевод битов в complex
	//string fileName;
	double Ps;
	complex<double> result;

public:

	/*signalSource(const string& fileName = "bits.txt") : fileName(fileName) {
		ifstream file(fileName, ios::in);
	}*/

	signalSource(const double& Ps = 2) : Ps(Ps) {
		result = { sqrt(Ps / 2.), 0.0 };
		srand(time(0));
	}

	complex<double> next() {
		return (rand() % 2 == 0 ? result : -result);
	}

};



class upsampling {
	// вставка нулей 2/4/8
	int n;
	int upsampling_n;
	complex<double> signalup;

public:

	upsampling(const int& upsampling_n = 1) : upsampling_n(upsampling_n) {
		n = 0;
	}

	complex<double> next(const complex<double>& signal) {
		if (n % upsampling_n == 0)
			signalup = signal;
		else
			signalup = 0;
		++n;
		return signalup;
	}

	int getUpsampling_n() {
		return upsampling_n;
	}

};



class FIR {
	// прием коэффициентов из файла
	// функция next для фильтрации
	double order;
	double coeff;
	complex<double> filteredSignal, temp;
	vector<double> coeffTable;
	deque<complex<double>> buff;
	string fileName;
	ifstream file;

public:

	FIR(const string& fileName = "filterTable.txt") : fileName(fileName), order(0), filteredSignal(0) {
		ifstream file(fileName, ios::in);

		if(file.is_open())
			while (!file.eof()) {
				file >> coeff;
				/*if (coeff == NULL)
					break;*/
				coeffTable.push_back(coeff);
				++order;
			}

		file.close();
		//buff.resize(order);
	}

	// цикл со сверткой с отдельной переменной которая отсчитывает свертку для коэфов
	// возвращает отфильтрованный отсчет
	complex<double> next(const complex<double>& signal) {
		if (buff.size() != order)
			buff.push_back(signal);
		else {
			buff.pop_front();
			buff.push_back(signal);
		}

		for (int i = 0; i < buff.size(); ++i) {
			filteredSignal += coeffTable[i] * buff[buff.size() - i - 1];
		}
		temp = filteredSignal;
		filteredSignal = {0.0, 0.0};

		return temp;
	}

	double getOrder() {
		return order;
	}

	void clear() {
		buff.clear();
	}

};



// измерение мощности сигнала
class PWR {
	double L, result, temp;
	deque<complex<double>> buff;

public:

	PWR(const double& L = 64) : L(L), result(0) {}

	// возвращает сумму квадратов амплитуд отсчетов, деленную на L
	double next(const complex<double>& signal) {
		if (buff.size() != L)
			buff.push_back(signal);
		else {
			buff.pop_front();
			buff.push_back(signal);
		}

		for (int i = 0; i < buff.size(); ++i) {
			result += pow(buff[i].real(), 2) + pow(buff[i].imag(), 2);
		}
		temp = result / L;
		result = 0;

		return temp;
	}

	double getPWR() {
		return temp;
	}


};



class awgn {
	// генерация шума next
	// pwr measure класс и сам шум
	default_random_engine generator;
	normal_distribution<double> distr;
	double Pn, upsampling_n;

public:

	awgn(const double& Pn = 2, const double& upsampling_n = 1) : Pn(Pn), upsampling_n(upsampling_n) {
		this->Pn *= upsampling_n;
		distr = normal_distribution<double>(0.0, 1.0);
	}

	double next(bool measured = 0, double snr = 0, double PWR = 0) {
		if (measured = 0)
			return std::sqrt(Pn / 2.) * distr(generator) / upsampling_n;
		else {
			Pn = pow(10., (snr - 10. * log10(PWR)) / 10.);
			return std::sqrt(Pn / 2.) * distr(generator) / upsampling_n;
		}
	}

};



class testAfterFilter {
	// убрать определенные биты так чтобы их не было в 0
	// т. о. проверить какие отсчеты являются информационными
};



class decision {
	// >= 0 = 1, < 0 = -1
	complex<double> result;
	double Ps;

public:

	decision(const double& Ps = 2) : Ps(Ps) {}

	// возвращает демодулированный отсчёт
	complex<double> next(const complex<double>& signal) {
		if (signal.real() >= 0.0)
			return { sqrt(Ps / 2.), 0.0 };
		else
			return {-sqrt(Ps / 2.), 0.0 };
	}

};



class biterror {
	// сравнение decision выхода и изначального signal через функцию next
	double N, correct;

public:

	biterror(const double& N = 1) : N(N), correct(0) {}

	// возвращает значение BER
	double next(complex<double> original, complex<double> final) {
		if (original == final)
			++correct;
		return (N - correct) / N;
	}

	void clear() {
		correct = 0;
	}

};




void writeComplexToPCMFile(const string& filename, complex<double> data) {
	ofstream file(filename, ios::app | ios::binary);
	if (file.is_open()) {
		int16_t real = static_cast<int16_t>(data.real() * INT16_MAX);
		int16_t imag = static_cast<int16_t>(data.imag() * INT16_MAX);
		file.write(reinterpret_cast<const char*>(&real), sizeof(int16_t));
		file.write(reinterpret_cast<const char*>(&imag), sizeof(int16_t));

		file.close();
		//cout << "Data written to PCM file successfully." << endl;
	}
	else {
		cout << "Unable to open file for writing." << endl;
	}
}




int main()
{
	double Ps = 2;
	double PsdB = 10. * log10(Ps);
	double Pn, PndB, BERRatio = 0, BERRatio_t = 0;

	int N = 10000;
	int upsampling_n = 2;
	int L = 64;

	complex<double> signal, signalProcessed;

	signalSource source(Ps);
	upsampling upsampler(upsampling_n);
	FIR shapingFilter("filterTable.txt");
	FIR matchedFilter("filterTable.txt");
	PWR measurer(L);
	decision decisionmaker(Ps);
	biterror biterr(N);
	biterror biterr_t(N);

	double delay = shapingFilter.getOrder();
	deque<complex<double>> buff_orig, buff_fin;

	ofstream out("out.txt");
	ofstream snrcheck("snr.txt");
	ofstream pwrout("pwr.txt");
	ofstream berout("ber.txt");

	for (int snr = 10; snr < 11; ++snr) {
		PndB = (PsdB - snr);
		Pn = pow(10., PndB / 10.);
		awgn noise(Pn, upsampling_n);
		for (int i = 0; i < N + delay; ++i) {
			// N + delay ?
			// if i == 0 запоминаем первый отсчет
			// if i >= delay находим
			// fifo размером delay, сравниваем первый и последний отсчеты
			signal = source.next();
			if (buff_orig.size() != delay + 1)
				buff_orig.push_back(signal);
			else {
				buff_orig.pop_front();
				buff_orig.push_back(signal);
			}
			for (int j = 0; j < upsampling_n; ++j) {
				signal = upsampler.next(signal);

				// вывод первоначального сигнала
				out << i << " " << j << "\t" << signal << "\t";

				// проверка мощности сигнала до шума
				double power = measurer.next(signal);

				signalProcessed = shapingFilter.next(signal);
				signalProcessed = { signalProcessed.real() + noise.next(1, snr, power), signalProcessed.imag() + noise.next(1, snr, power) };

				// вывод для проверки осш в audition
				snrcheck << signalProcessed.real() << endl;
				//writeComplexToPCMFile("output.pcm", signalProcessed);
				//writeComplexToPCMFile("output.pcm", noise.next());

				signalProcessed = matchedFilter.next(signalProcessed);

				// вывод для проверки мощности сигнала после шума
				//pwrout << signalProcessed << "\t" << measurer.next(signalProcessed) << endl;

				// вывод для проверки на выбор нужного бита

				signalProcessed = decisionmaker.next(signalProcessed);

				// вывод конечного сигнала для проверки delay : соотвествует order фильтра
				out << signalProcessed << endl;
				BERRatio = biterr.next(signal, signalProcessed);
			}

			/*if (buff_fin.size() != 1)
				buff_fin.push_back(signalProcessed);
			else {
				buff_fin.pop_front();
				buff_fin.push_back(signalProcessed);
			}*/

			// конечный бит после апсемплинга должен нести информацию
			if(i >= delay)
				BERRatio_t = biterr_t.next(buff_orig[delay], signalProcessed);
		}

		berout << snr << "\t" << BERRatio << "\t" << BERRatio_t << "\t" << 0.5 * erfc(sqrt(pow(10., snr / 10.))) << endl;
		
		shapingFilter.clear();
		matchedFilter.clear();
		biterr.clear();
	}


	out.close();
	snrcheck.close();
	pwrout.close();
	berout.close();

	return 0;
}
