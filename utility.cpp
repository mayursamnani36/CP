//order_by_key gives the lower bound of the given element find_by_order gives iterator to
//the position of the element so * it to get the position. If elements are not distinct
//use pair<int, int> instead of int and make second value unique
#include<ext/pb_ds/assoc_container.hpp>
#include<ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> PBDS;


lli kfib(lli n) {return (pow((1 + sqrt(5)), n) - pow((1 - sqrt(5)), n)) / (sqrt(5) * pow(2, n));}

void sieve() {
	// O(nlog(log(n)))
	lli maxN = 1000000;
	for (lli i = 3 ; i <= maxN ; i += 2) {primeArr[i] = true;}
	primeArr[2] = true;
	primes.pb(2);
	for (lli i = 3 ; i <= maxN ; i += 2) {
		if (primeArr[i]) {
			primes.pb(i);
			for (lli j = i * i ; j <= maxN ; j += i) {primeArr[j] = false;}
		}
	}
}

vector<bool> s_sieve(lli m, lli n) {
	vector<bool> segment(n - m + 1, true);
	for (auto x : primes) {
		if (x * x > n) {break;}
		lli start = (m / x) * x;
		if (x >= m &&  x <= n) {start = 2 * x;}
		for (lli i = start; i <= n ; i += x) {
			if (i - m >= 0) segment[i - m] = false;
		}
	}
	// for (i = m ; i <= n ; i++) {
	// 		if (segment[i - m] && i != 1) {cout << i << nl;}
	// }
	return segment;
}

lli ETF(lli n) {
	// Euler Totient Function (no. of no.s from 1 to n which are coprime to n) O(sqrt(n))
	lli res = n;
	for (lli i = 2 ; i * i <= n ; i++) {
		if (n % i == 0) {
			res /= i;
			res *= (i - 1);
			while (n % i == 0) {n /= i;}
		}
	}
	if (n > 1) {
		res /= n;
		res *= (n - 1);
	}
	return res;
}

lli binpowmod(lli a, lli b, lli m) {
	a %= m;
	if (b == 0) {return 1;}
	lli res = 1;
	while (b > 0) {
		if (b & 1) {res = (res * a) % m;}
		b >>= 1;
		a = a * a % m;
	}
	return res;
}

lli binpow(lli a, lli b) {
	// O(log(n))
	lli ans = 1;
	while (b > 0) {
		if (b & 1) {ans *= a;}
		b >>= 1;
		a *= a;
	}
	return ans;
}

bool isPrime(lli num) {
	// O(sqrt(n))
	bool flag = true;
	for (lli i = 2; i <= sqrt(num); i++) {
		if (num % i == 0) {
			flag = false;
			break;
		}
	}
	return flag;
}

struct cmp {
	bool operator() (const pair<lli, lli> &a, const pair<lli, lli> &b) const {
		lli lena = a.second - a.first + 1;
		lli lenb = b.second - b.first + 1;
		if (lena == lenb) return a.first < b.first;
		return lena > lenb;
	}
};

string baseConvertor(lli sb, lli db, lli num) {
	lli bt = 0;
	lli p = 0;
	while (num) {
		lli x = num % 10;
		num /= 10;
		bt += binpow(sb, p) * x;
		p++;
	}
	string ans = "";
	while (bt) {
		lli x = bt % db;
		bt /= db;
		ans = to_string(x) + ans;
	}
	return ans;
}

int csb(int n) {
	int count = 0;
	while (n) {
		count++;
		n = n & (n - 1); // remove set bits from right to left
		//__builtin_popcount(n) also gives the answer
	}
	return count;
}

void extendedEuclid(lli a, lli b) {
	//solx is the modulo inverse of a wrt to b if gcd(a,b)==1
	//Also it can be negative so make sure to add b if needed

	lli solx, soly, GCD; // copy these outside
	if (b == 0) {
		solx = 1;
		soly = 0;
		GCD = a;
		return;
	}
	extendedEuclid(b, a % b);
	lli cx = soly;
	lli cy = solx - (a / b) * soly;
	solx = cx;
	soly = cy;
}

void modfact(int n) {
	fact[0] = 1;
	for (int i = 1 ; i <= n; i++)
		fact[i] = (fact[i - 1] * i) % mod;
}

lli C(lli n, lli r) {
	if (n < r) return 0;
	if (n == r || r == 0) return 1;
	if (r == 1 || r == n - 1) return n % mod;
	extendedEuclid(fact[r], mod);
	lli temp1 = (solx + mod) % mod;
	extendedEuclid(fact[n - r] % mod, mod);
	lli temp2 = (solx + mod) % mod;
	return ((fact[n] * temp1) % mod * temp2 % mod) % mod;
}


// For matrix expo
vector<vector<lli> > multiply(vector<vector<lli>> A , vector<vector<lli>> B) {
	vector<vector<lli>> C(n, vector<lli>(n));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int x = 0; x < n; x++) {
				C[i][j] = (C[i][j] + (A[i][x] * B[x][j]) % MOD) % MOD;
			}
		}
	}
	return C;
}
vector<vector<lli> >  pow(vector<vector<lli> > A, lli p) {
	if (p == 1) {return A;}
	if (p & 1) {return multiply(A, pow(A, p - 1));}
	else {
		vector<vector<lli> > X = pow(A, p / 2);
		return multiply(X, X);
	}
}
//matrix expo end