//============================================================================
// Name        : Module_2.cpp
// Author      : wufy
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <thread>
#include <mutex>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
using namespace std;

#ifdef _WIN32
#include <windows.h>

void sleep(unsigned milliseconds)
{
	Sleep(milliseconds);
}
#else
#include <unistd.h>

void sleep(unsigned milliseconds)
{
	usleep(milliseconds * 1000); // takes microseconds
}
#endif


#define CAPACITY 12303.00
#define ALLCARNUM 15752.00
#define ALLCARNUM_INT 15752
#define Vkf 4.0
#define Vs  8.0

#define Vkf_W 0.5
#define Vs_W  1.0

#define WALK_CAPACITY 3717.00

const double eps = 1e-7;


// 记录NormalStats已经算过的值，提高效率
std::map<std::string, double> mapVis;

class People
{
public:
	People(double UseCar = 0.80, int Threshold = 10)
	{
		UseCar_ = UseCar;
		Threshold_ = Threshold;
		carTimes_ = 0;
		walkTimes_ = 0;

		vauleCarVec.clear();
		vauleWalkVec.clear();

		valueCarExcept = 0.00;
		valueCarVariance = 0.00;

		valueWalkExcept = 0.00;
		valueWalkVariance = 0.00;

	}
	void correct(double traffic)
	{
		// 效用函数

		//		double Uj = getUj(traffic);
		//		double Us = -0.154 / 1.5 + 0.612  * 4.79;
		//
		//		double Uj_e = exp(Uj);
		//		double Us_e = exp(Us);
		//
		//		double UseCar_new =  Uj_e / (Uj_e + Us_e);
		//
		//		double step = getSteps(traffic, ALLCARNUM * UseCar_new);
		//
		//		// cout << "step:" << step << endl;
		//
		//		UseCar_ =  UseCar_  * (1 - step)  + step * UseCar_new ;

		if (isUserCar_)
		{
			double Uj = getUj(traffic);
			vauleCarVec.push_back(Uj);
			valueCarExcept = (valueCarExcept * (carTimes_ - 1) + Uj) / carTimes_;

			double ret = 0;
			for (size_t i = 0; i < vauleCarVec.size(); i++)
			{
				ret += (vauleCarVec[i] - valueCarExcept) * (vauleCarVec[i] - valueCarExcept);
			}
			valueCarVariance = ret / vauleCarVec.size();
		}
		else
		{
			double Us = getUs(ALLCARNUM - traffic);
			vauleWalkVec.push_back(Us);
			valueWalkExcept = (valueWalkExcept * (walkTimes_ - 1) + Us) / walkTimes_;

			double ret = 0;
			for (size_t i = 0; i < vauleWalkVec.size(); i++)
			{
				ret += (vauleWalkVec[i] - valueWalkExcept) * (vauleWalkVec[i] - valueWalkExcept);
			}
			valueWalkVariance = ret / vauleWalkVec.size();
		}


		valueCarVariance = max(0.1, valueCarVariance);
		valueWalkVariance = max(0.1, valueWalkVariance);

		if (valueCarVariance > eps && valueWalkVariance > eps && valueCarExcept > eps && valueWalkExcept > eps)
		{

			double minVCar = valueCarExcept - 3.00 * sqrt(valueCarVariance);
			double maxVCar = valueCarExcept + 3.00 * sqrt(valueCarVariance);
			double minVWalk = valueWalkExcept - 3.00 * sqrt(valueWalkVariance);
			double maxVWalk = valueWalkExcept + 3.00 * sqrt(valueWalkVariance);

			double P1Car = NormalStats(minVCar, Vkf, valueCarExcept, valueCarVariance);
			double P2Car = NormalStats(Vkf, Vs, valueCarExcept, valueCarVariance);
			double P3Car = NormalStats(Vs, maxVCar, valueCarExcept, valueCarVariance);

			double P1Walk = NormalStats(minVWalk, Vkf_W, valueWalkExcept, valueWalkVariance);
			double P2Walk = NormalStats(Vkf_W, Vs_W, valueWalkExcept, valueWalkVariance);
			double P3Walk = NormalStats(Vs_W, maxVWalk, valueWalkExcept, valueWalkVariance);

			double Vk1Car = 0.00;
			double Vk2Car = 0.00;
			double Vk3Car = 0.00;
			if (P1Car > 0.00) {
				Vk1Car = getValueWithProbability(minVCar, Vkf, valueCarExcept, valueCarVariance) / P1Car;
			}
			if (P2Car > 0.00) {
				Vk2Car = getValueWithProbability(Vkf, Vs, valueCarExcept, valueCarVariance) / P2Car;
			}
			if (P3Car > 0.00) {
				Vk3Car = getValueWithProbability(Vs, maxVCar, valueCarExcept, valueCarVariance);
			}

			double Vk1Walk = 0.00;
			double Vk2Walk = 0.00;
			double Vk3Walk = 0.00;

			if (P1Walk > 0.00) {
				Vk1Walk = getValueWithProbability(minVWalk, Vkf_W, valueWalkExcept, valueWalkVariance) / P1Walk;
			}
			if (P2Walk > 0.00) {
				Vk2Walk = getValueWithProbability(Vkf_W, Vs_W, valueWalkExcept, valueWalkVariance) / P2Walk;
			}
			if (P3Walk > 0.00) {
				Vk3Walk = getValueWithProbability(Vs_W, maxVWalk, valueWalkExcept, valueWalkVariance) / P3Walk;
			}

			double PECar1 = Vk1Car * getWeights(P1Car, 0.69);
			double PECar2 = Vk2Car * getWeights(P2Car, 0.61);
			double PECar3 = Vk3Car * getWeights(P3Car, 0.69);

			double PEWalk1 = Vk1Walk * getWeights(P1Walk, 0.69);
			double PEWalk2 = Vk2Walk * getWeights(P2Walk, 0.61);
			double PEWalk3 = Vk3Walk * getWeights(P3Walk, 0.69);

			double PECar = PECar1 + PECar2 + PECar3;
			double PEWalk = PEWalk1 + PEWalk2 + PEWalk3;

			UseCar_ = exp(PECar) / (exp(PECar) + exp(PEWalk));
			// UseCar_ = PECar / (PECar + PEWalk);

		}
	}

	bool IsUseCar()
	{
		double randD = 0.01 * (rand() % 100);
		isUserCar_ = randD < UseCar_;
		if (isUserCar_)
			carTimes_++;
		else
			walkTimes_++;

		return isUserCar_;
	}

	double UseCar_;

private:

	int carTimes_;
	int walkTimes_;
	int Threshold_;
	bool isUserCar_;

	double valueCarExcept;
	double valueCarVariance;

	double valueWalkExcept;
	double valueWalkVariance;

	vector<double> vauleCarVec;
	vector<double> vauleWalkVec;

	double getTime(double traffic)
	{
		return (1 + 0.15 * (traffic / CAPACITY) * (traffic / CAPACITY)  * (traffic / CAPACITY)  * (traffic / CAPACITY)) / 11.00;
	}

	//	double getSteps(double traffic, double traffic_new)
	//	{
	//		// return 1.00 / (double)times;
	//
	//		//		double ans = 0.00;
	//		//		double ret = 10000;
	//		//
	//		//		double temp_ans = 0.00;
	//		//		for (int i = 0; i < 100; i++)
	//		//		{
	//		//			// cout << "Uj: " <<  getUj(traffic + temp_ans * (traffic_new - traffic)) << endl;
	//		//			// cout << "Us: " << getUs((ALLCARNUM - traffic) + temp_ans * (traffic - traffic_new)) << endl;
	//		//			double temp = (traffic_new - traffic) * getUj(traffic + temp_ans * (traffic_new - traffic));
	//		//			temp += (traffic - traffic_new) * getUs((ALLCARNUM - traffic) + temp_ans * (traffic - traffic_new));
	//		//
	//		//			// cout << "temp_step" << i << ": " << temp << endl;
	//		//
	//		//			if (abs(temp) < ret)
	//		//			{
	//		//				ans = temp_ans;
	//		//				ret = abs(temp);
	//		//			}
	//		//			temp_ans += 0.01;
	//		//		}
	//		//
	//		//		return ans;
	//	}

	double getUj(double traffic)
	{
		double Tj = getTime(traffic);
		double Fj = 0.20;
		double Cj = CAPACITY / traffic*2.67;
		double Fj1 = (2 * (Cj - Tj) - (Cj - Fj)) / (Cj - Tj);

		return -0.154 - 0.202 * Fj1 + 0.612 *Cj;
	}

	double getUs(double traffic)
	{
		return -0.154 / 1.5 + 0.612  *  WALK_CAPACITY / traffic * 4.67;
	}

	double getWeights(double P, double r)
	{
		return pow(P, r) / pow(pow(P, r) + pow(1 - P, r), 1 / r);
	}


	double getValue(double Vkn)
	{
		if (Vkn < Vkf)
		{
			return -2.25 * pow((Vkf - Vkn), 0.88);
		}
		else if (Vkn <= Vs)
		{
			return pow((Vkn - Vkf), 0.88);
		}
		else
		{
			return pow((Vs - Vkf), 0.88);
		}
	}


	double NormalExp(double x, double expect = 0.0, double variance = 1.0, bool withVkn = false)
	{
		if (withVkn)
		{
			return exp(-(x - expect) * (x - expect) / (2 * variance)) * getValue(x);
		}
		else
		{
			return exp(-(x - expect) * (x - expect) / (2 * variance));
		}
	}

	double NormalStats(double a, double b, double expect = 0.0, double variance = 1.0, bool withVkn = false)
	{
		if (a >= b)
			return 0.00;

		double s1 = 0;
		double s2 = (b - a) * (NormalExp(a, expect, variance, withVkn) + NormalExp(b, expect, variance, withVkn)) / 2;

		if (s2 < eps)
			return 0.00;

		std::string strKey = std::to_string(a) + std::to_string(b) + std::to_string(expect) + std::to_string(variance) + std::to_string(withVkn);

		std::map<std::string , double > ::iterator iter = mapVis.find(strKey);
		if (iter != mapVis.end())
			return iter->second;


		for (int i = 1; fabs(s1 - s2) > eps; i *= 2)
		{
			double h = (b - a) / (double)i;
			s1 = s2;
			s2 = 0;
			for (int j = 0; j < i; j++)
			{
				s2 += h * NormalExp(a + ((double)j + 0.5) * h, expect, variance, withVkn);
			}
			s2 = (s1 + s2) / 2;
		}
		int ans = s2 * sqrt(1 / (8 * atan(1.0) * variance));
	    mapVis.insert(make_pair(strKey, ans));
		return ans;
	}


	double getValueWithProbability(double a, double b, double expect, double variance)
	{

		if (a >= b)
			return 0.00;

		double steps = (b - a) / 100;

		double sum = 0.0;
		double start = a;
		double end = a + steps;

		for (int i = 0; i < 100; i++)
		{
			sum += NormalStats(start, end, expect, variance) * (getValue(start) + getValue(end)) / 2.00;
			start = end;
			end += steps;
		}

		return sum;
	}


};

static People people[ALLCARNUM_INT + 9];
static int finishedThreadNum = 0;
static std::mutex finishMutex;

void correctThread(int id, int traffic)
{
	for (int i = id; i < ALLCARNUM_INT; i += 10)
	{
		people[i].correct(traffic);
	}
	{
		std::lock_guard<std::mutex> lck(finishMutex);
		finishedThreadNum++;
	}
}




int main()
{
	// cout << NormalStats(-1 , 1, 2 , 2) << endl;

	for (int i = 0; i < ALLCARNUM_INT; i++)
	{
		people[i] = People();
	}

	int day = 1;
	while (1)
	{
		int traffic = 0;
		for (int i = 0; i < ALLCARNUM_INT; i++)
		{
			if (people[i].IsUseCar())
				traffic++;
		}

		traffic = (double)traffic / 2.80;

		// traffic = people[0].UseCar_ * ALLCARNUM;
		for (int i = 0; i < 10; i++)
		{
			std::thread thr(&correctThread, i, traffic);//创建一个分支线程，回调到myThread函数里
			thr.detach();
		}
		//for (int i = 0; i < 15752; i++)
		//{
			// cout << "day: " << day << "people: " << i <<endl;
			//people[i].correct(traffic);
		//}
		// cout << "第" << day++ << "天的交通量是" << traffic << endl;
		while (1)
		{
			{
				std::lock_guard<std::mutex> lck(finishMutex);
				if (finishedThreadNum >= 10)
					break;
			}
			sleep(100);
		}
		cout << "day: " << day << "  traffic:" << traffic << endl;

		day++;
		// getchar();
		if (day > 200)
			break;

	}


	return 0;
}
