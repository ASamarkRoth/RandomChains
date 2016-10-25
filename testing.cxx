#include <iostream>
//#include "TH1.h"

using namespace std;

template <unsigned const int hej, unsigned const int da>
class rectangle {
	private:
	int length, width;
	int* hejda;

	public:
	rectangle(int w, int l);
	int get_width();
	int get_length();
	void set_width(int w);
	void set_length(int l);

	int get_area();


};

template <unsigned const int hej, unsigned const int da>
rectangle<hej,da>::rectangle(int w, int l) {
	width = w;
	length = l;
	hejda = new int[hej][da];
}

template <unsigned const int hej, unsigned const int da>
int rectangle<hej,da>::get_width() {
	return width;
}
template <unsigned const int hej, unsigned const int da>
int rectangle<hej,da>::get_length() {
	return length;
}

template <unsigned const int hej, unsigned const int da>
void rectangle<hej,da>::set_width(int w){
	width = w;
}
template <unsigned const int hej, unsigned const int da>
void rectangle<hej,da>::set_length(int l){
	length = l;
}
template <unsigned const int hej, unsigned const int da>
int rectangle<hej,da>::get_area() {
	return width*length;
}

int a_main() {
	
	rectangle<5,5> rec = rectangle<5,5>(10, 10);

	cout << "This is width " << rec.get_width() << " and this is length " << rec.get_length() << endl;

	/*
	array<double, 2> hihi;
	hihi[0] = 3;
	hihi[1] = 2;
	
	
	for(int j = 0; j < 2; j++) {
		cout << "hihi = " << hihi[j] << endl;
	}
	*/
	
	return 0;
	}
