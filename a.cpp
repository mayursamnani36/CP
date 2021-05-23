#pragma GCC optimize("trapv")
#pragma GCC optimize("Ofast")
#pragma GCC optimization ("unroll-loops")
#pragma GCC target ("avx2")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")
#include "bits/stdc++.h"
using namespace std;
#ifdef ONLINE_JUDGE
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> OST;
#endif
#define int long long
#define double long double
#define mp make_pair
#define fi first
#define se second
#define pb push_back
#define all(a) a.begin(),a.end()
const int MOD = 1000000007;
const int MOD2 = 998244353;
const int BIG = 1197423052705914509LL;
mt19937 rng((int) std::chrono::steady_clock::now().time_since_epoch().count());
const int MAXN = 4e5 + 10;
struct custom_hash {
	static uint64_t splitmix64(uint64_t x) {
		x += 0x9e3779b97f4a7c15;
		x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
		x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
		return x ^ (x >> 31);
	}
	size_t operator()(uint64_t a) const {
		static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
		return splitmix64(a + FIXED_RANDOM);
	}
	template<class T> size_t operator()(T a) const {
		static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
		hash<T> x;
		return splitmix64(x(a) + FIXED_RANDOM);
	}
	template<class T, class H> size_t operator()(pair<T, H> a) const {
		static const uint64_t FIXED_RANDOM = chrono::steady_clock::now().time_since_epoch().count();
		hash<T> x;
		hash<H> y;
		return splitmix64(x(a.f) * 37 + y(a.s) + FIXED_RANDOM);
	}
};
template<class T, class H>using umap = unordered_map<T, H, custom_hash>;
void debug(int x) {}
void solve(int test_case) {
	//Checker solution for small cases
	int n = rand() % 9 + 2;
	int a[n + 1];
	for (int i = 1; i <= n; i++)a[i] = rand() % 3 + 1;
	bool b[n + 1];
	for (int i = 1; i <= n; i++)b[i] = 0;
	int dpa = 0;
	for (int i = 0; i <= n; i++) {
		if (i) {
			sort(b + 1, b + n + 1);
			b[1] = 1;
			sort(b + 1, b + n + 1);
		}
		do {
			vector<int>v0, v1;
			for (int j = 1; j <= n; j++) {
				if (b[j] == 0)v0.pb(a[j]);
				else v1.pb(a[j]);
			}
			int cnt = 0;
			for (int j = 0; j < v0.size(); j++) {
				if (j == 0 || v0[j] != v0[j - 1])cnt++;
			}
			for (int j = 0; j < v1.size(); j++) {
				if (j == 0 || v1[j] != v1[j - 1])cnt++;
			}
			dpa = max(dpa, cnt);

		} while (next_permutation(b + 1, b + n + 1));
	}
	//Insert target participant solution

	if (dpa != ans) {
		cout << "WA on test case " << test_case << "\n";
		cout << n << "\n";
		for (int i = 1; i <= n; i++)cout << a[i] << " ";
		cout << "\n";
		cout << "Expected " << dpa << ", found " << ans << "\n";
	}

}
signed main() {
	time_t t = clock();
	srand(time(NULL));
	ios::sync_with_stdio(0);
	cin.tie(0);
	int T = 100; //cin >> T;
	int test_case = 1;
	while (T--) {
		solve(test_case);
		test_case++;
	}
	cerr << "Program terminated successfully\n";
	t = clock() - t;
	cerr << "Time used: " << fixed << setprecision(3) << (double)(t * 1.0 / CLOCKS_PER_SEC) << " seconds\n";
}