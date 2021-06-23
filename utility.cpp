# Stay calm
# Think reverse
# Think DP
# Think Graph
# Check i++ i--
# Read question again
# Think digit by digit

//order_by_key gives the lower bound of the given element find_by_order gives iterator to
//the position of the element so * it to get the position. If elements are not distinct
//use pair<int, int> instead of int and make second value unique
#include<ext/pb_ds/assoc_container.hpp>
#include<ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> ordered_set;

lli getInvCount(vector<lli> &arr) {
	lli key;
	ordered_set set1;
	set1.insert(arr[0]);
	lli invcount = 0;
	lli n = arr.size();
	for (lli i = 1; i < n; i++) {
		set1.insert(arr[i]);
		key = set1.order_of_key(arr[i] + 1);
		invcount += set1.size() - key;
	}
	return invcount;
}

// Numbers whose mask is same when multiplied give perfect square
lli mask(lli n) {
	lli ans = 1;
	lli count = 0;
	while (n % 2 == 0) {
		count++;
		n /= 2;
	}
	if (count & 1) {ans *= 2;}
	for (lli i = 3; i * i <= n; i += 2) {
		count = 0;
		while (n % i == 0) {
			count++;
			n = n / i;
		}
		if (count & 1) {ans *= i;}
	}
	if (n > 1) {ans *= n;}
	return ans;
}

//srand(atoi(argv[1])); put this in main as well
lli rand(lli a, lli b) {return a + rand() % (b - a + 1);}

lli findPivot(vector<lli> &arr) {
	lli s = 0;
	lli e = arr.size() - 1;
	if (arr[s] < arr[e]) {return e;}
	lli m;
	while (s <= e) {

		m = (s + e) / 2;
		if ((m < e && arr[m] > arr[m + 1]) || (m > s && arr[m] < arr[m - 1])) {return m;}
		if (arr[s] >= arr[m]) {e = m - 1;}
		else {s = m + 1;}
	}
	return m;
}

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
	// for (lli i = m ; i <= n ; i++) {
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

lli CRT(vector<lli> mods, vector<lli> rems) {
	lli k = mods.size();
	lli p = 1;
	lli i;
	for (i = 0 ; i < k ; i++) {p *= mods[i];}
	vector<lli> pp(k), inv(k);

	for (i = 0 ; i < k ; i++) {
		pp[i] = p / mods[i];
		extendedEuclid(pp[i], mods[i]);
		inv[i] = (solx + mods[i]) % mods[i];
	}

	lli ans = 0;
	for (i = 0 ; i < k ; i++) {
		ans += ((rems[i] * pp[i] * inv[i]) % p);
	}
	return ans % p;
}

void modfact(lli n) {
	fact[0] = 1;
	for (lli i = 1 ; i <= n; i++)
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

// binomial coeff
vector<vector<lli>> nCr(61, vector<lli>(61, 0));
for (int i = 0; i < 61; i++) {
	for (int j = 0; j <= i; j++) {
		if (i == j || j == 0) {nCr[i][j] = 1;}
		else {nCr[i][j] = nCr[i - 1][j - 1] + nCr[i - 1][j];}
	}
}

// For matrix expo---------------------------------------------------------------------------
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
vector<vector<lli> > pow(vector<vector<lli> > A, lli p) {
	if (p == 1) {return A;}
	if (p & 1) {return multiply(A, pow(A, p - 1));}
	else {
		vector<vector<lli> > X = pow(A, p / 2);
		return multiply(X, X);
	}
}
//matrix expo end----------------------------------------------------------------------------


// Segment Tree implementation for min range query and update ----------------------------------------------------
// complete overlap matlab node ki range lies in query ki range

void buildTree(vector<lli> &v, lli ss, lli se, vector<lli>&tree, lli index) {
	// at a particular index we store min from ss to se
	if (ss == se) {tree[index] = v[ss]; return;}
	lli mid = ss + (se - ss) / 2;
	buildTree(v, ss, mid, tree, 2 * index);
	buildTree(v, mid + 1, se, tree, 2 * index + 1);
	tree[index] = min(tree[2 * index], tree[2 * index + 1]);
	return;
}

lli query(vector<lli> &tree, lli ss, lli se, lli qs, lli qe, lli index) {
	// Complete overlap
	if (ss >= qs && se <= qe) {return tree[index];}
	// No overlap
	if (qe < ss || qs > se) {return INT_MAX;}
	// Partial overlap
	lli mid = ss + (se - ss) / 2;
	lli left = query(tree, ss, mid, qs, qe, 2 * index);
	lli right = query(tree, mid + 1, se, qs, qe, 2 * index + 1);
	return min(left, right);
}

void updateNode(vector<lli> &tree, lli ss, lli se, lli i , lli increment, lli index) {
	// Case i is out of bounds
	if (i < ss || i > se) {return;}
	if (ss == se) {tree[index] += increment; return;}
	lli mid = ss + (se - ss) / 2;
	updateNode(tree, ss, mid, i, increment, 2 * index);
	updateNode(tree, mid + 1, se, i, increment, 2 * index + 1);
	tree[index] = min(tree[2 * index], tree[2 * index + 1]);
	return;
}

void updateRange(vector<lli> &tree, lli ss, lli se, lli l, lli r, lli increment, lli index) {
	//Out of bounds
	if (l > se || r < ss) {return;}
	//leaf node (dont put it before out of bounds)
	if (ss == se) {
		tree[index] += increment;
		return;
	}
	lli mid = ss + (se - ss) / 2;
	updateRange(tree, ss, mid, l, r, increment, 2 * index);
	updateRange(tree, mid + 1, se, l, r, increment, 2 * index + 1);
	tree[index] = min(tree[2 * index], tree[2 * index + 1]);
	return;
}

void updateRangeLazy(vector<lli> &tree, lli ss, lli se, lli l, lli r, lli increment, lli index) {
	//before going down resolve the value if exists
	if (lazy[index] != 0) {
		tree[index] += lazy[index];
		//non leaf node
		if (ss != se) {
			lazy[2 * index] += lazy[index];
			lazy[2 * index + 1] += lazy[index];
		}
		lazy[index] = 0;
	}

	//no overlap
	if (l > se || r < ss) {return;}

	//complete overlap
	if (ss >= l && se <= r) {
		tree[index] += increment;
		if (ss != se) {
			lazy[2 * index] += increment;
			lazy[2 * index + 1] += increment;
		}
		return;
	}

	lli mid = ss + (se - ss) / 2;
	updateRangeLazy(tree, ss, mid, l, r, increment, 2 * index);
	updateRangeLazy(tree, mid + 1, se, l, r, increment, 2 * index + 1);
	tree[index] = min(tree[2 * index], tree[2 * index + 1]);
	return;
}

lli queryLazy(vector<lli> &tree, lli ss, lli se, lli qs, lli qe, lli index) {
	// resolve the lazy value at current node
	if (lazy[index] != 0) {
		tree[index] += lazy[index];
		//non leaf node
		if (ss != se) {
			lazy[2 * index] += lazy[index];
			lazy[2 * index + 1] += lazy[index];
		}
		lazy[index] = 0;
	}

	// No overlap
	if (ss > qe || se < qs) {return INT_MAX;} // for range mininum query

	// Complete overlap
	if (ss >= qs && se <= qe) {return tree[index];}

	lli mid = ss + (se - ss) / 2;
	lli left = queryLazy(tree, ss, mid, qs, qe , 2 * index);
	lli right = queryLazy(tree, mid + 1, se, qs, qe, 2 * index + 1);
	return min(left, right);
}


// Segment tree end ---------------------------------------------------------------------


// Fenwick tree --------------------------------------------------------------------------

void update(lli i, lli inc, lli n) {
	while (i <= n) {
		BIT[i] += inc;
		i += (i & (-i));
	}
}

lli query(lli i) {
	lli sum = 0;
	while (i > 0) {
		sum += BIT[i];
		i -= (i & (-i));
	}
	return sum;
}

// Fenwick tree end --------------------------------------------------------------------------


//DSU-----------------------------------------------------------------------------------------

// make a parent array with each element as parent of its own i.e parent[i] = i;
// initialize size array as 1
fo(i, n + 1) {parent[i] = i; sz[i] = 1;}

lli parent[1000001];
lli sz[1000001];

lli find_set(lli v) {
	if (v == parent[v]) {return v;}
	return parent[v] = find_set(parent[v]);
}

void union_sets(lli a, lli b) {
	a = find_set(a);
	b = find_set(b);
	if (a != b) {
		if (sz[a] < sz[b]) {swap(a, b);}
		parent[b] = a;
		sz[a] += sz[b];
	}
}

// for finding different connected componenets
set<lli> st;
for (i = 1 ; i <= n ; i++) {
	lli temp = i;
	while (parent[temp] != temp) {
		temp = parent[temp];
	}
	st.insert(parent[temp]);
}


//DSU end ------------------------------------------------------------------------------------

void dijkstra(vector<vector<plli>> &adj, vector<lli> &dist, vector<bool> &vis, lli start) {
	// adj is node -> node, dis
	// st is dist, node
	set<plli> st;
	dist[start] = 0;
	st.insert({0, start});
	while (st.size()) {
		plli p = *st.begin();
		st.erase(st.begin());
		lli dis = p.ff;
		lli node = p.ss;
		vis[node] = true;
		for (auto nbr : adj[node]) {
			lli target = nbr.ff;
			lli tempdist = nbr.ss;
			if (dis + tempdist < dist[target]) {
				auto f = st.find({dist[target], target});
				if (f != st.end()) {
					st.erase(f);
				}
				dist[target] = dis + tempdist;
			}
			if (!vis[target])st.insert({dist[target], target});
		}
	}
}

void floydWarshall(vector<vector<lli>> &adjMat) {
	lli n = adjMat.size();
	lli i, j, k;
	//Phases, in kth phase we included vertices (1, 2 ... k) as intermediate vertex
	for (k = 0 ; k < n ; k++) {
		//Iterate over entire matrix
		for (i = 0 ; i < n ; i++) {
			for (j = 0 ; j < n ; j++) {
				adjMat[i][j] = min(adjMat[i][j], adjMat[i][k] + adjMat[k][j]);
			}
		}
	}
}

vector<lli> bellmanFord(lli n, lli src, vector<vector<lli>> &edgeList) {
	vector<lli> dist(n + 1, INT_MAX);
	dist[src] = 0;

	// relax all edges n-1 times
	for (lli i = 0 ; i < n - 1 ; i++) {
		for (auto edge : edgeList) {
			lli u = edge[0];
			lli v = edge[1];
			lli wt = edge[2];
			if (dist[u] != INT_MAX && dist[u] + wt < dist[v]) {dist[v] = dist[u] + wt;}
		}
	}

	// for -ve wt cycle
	for (auto edge : edgeList) {
		lli u = edge[0];
		lli v = edge[1];
		lli wt = edge[2];
		if (dist[u] != INT_MAX && dist[u] + wt < dist[v]) {return {};}
	}

	return dist;
}