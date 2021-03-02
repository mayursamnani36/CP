#include<bits/stdc++.h>
using namespace std;

int rand(int a, int b) {return a + rand() % (b - a + 1);}

int main(int argc, char* argv[]) {
	srand(atoi(argv[1]));
	int n = rand(1, 1000);
	cout << n << nl;
	for (i = 0 ; i <= n ; i++) {
		int x;
		x = rand(1, n);
		cout << x << " ";
	}
	cout << nl;
	return 0;
}