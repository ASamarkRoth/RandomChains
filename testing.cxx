#include <iostream>
//#include "TH1.h"

using namespace std;

class rectangle {
	private:
	int length, width;

	public:
	rectangle(int w, int l);
	int get_width();
	int get_length();
	void set_width(int w);
	void set_length(int l);

	int get_area();


};

rectangle::rectangle(int w, int l) {
	width = w;
	length = l;
}

int rectangle::get_width() {
	return width;
}
int rectangle::get_length() {
	return length;
}

void rectangle::set_width(int w){
	width = w;
}
void rectangle::set_length(int l){
	length = l;
}
int rectangle::get_area() {
	return width*length;
}

int a_main() {
	
	rectangle rec = rectangle(10, 10);

	cout << "This is width " << rec.get_width() << " and this is length " << rec.get_length() << endl;

	array<double, 2> hihi;
	hihi[0] = 3;
	hihi[1] = 2;
	
	for(int j = 0; j < 2; j++) {
		cout << "hihi = " << hihi[j] << endl;
	}
	
	return 0;
	}
