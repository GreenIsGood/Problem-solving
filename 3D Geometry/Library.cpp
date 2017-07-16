#define _USE_MATH_DEFINES
#include "bits/stdc++.h"

using namespace std;

typedef complex<double> point;
#define X real()
#define Y imag()
#define vec(a,b) ((b)-(a))
#define dot(a,b) ((conj(a)*(b)).real())
#define lengthSqr(v) (dot(v,v))
#define same(a,b) (lengthSqr(vec(a,b))<EPS)
#define EPS 1e-9

// Distance between two points on sphere
double dist(double lat1, double lng1, double lat2,
		double lng2, double R) {
	double Dlat = lat2 - lat1;
	double Dlng = lng2 - lng1;
	double A = sin(Dlat / 2) * sin(Dlat / 2)
			+ cos(lat1) * cos(lat2) * sin(Dlng / 2) * sin(Dlng / 2);
	double C = 2 * atan2(sqrt(A), sqrt(1 - A));
	return R * C;
}

double ONE = 1;
struct point3D {
	double v[3];
	point3D() {
		for (int i = 0; i < 3; ++i) {
			this->v[i] = 0;
		}
	}
	point3D(const double v[3]) {
		for (int i = 0; i < 3; ++i) {
			this->v[i] = v[i];
		}
	}
	point3D(double lon, double lat, double r) {
		x() = r * cos(lat) * cos(lon);
		y() = r * cos(lat) * sin(lon);
		z() = r * sin(lat);
	}

	double& operator [](int idx) {
		return idx < 3 ? v[idx] : (ONE = 1);
	}
	double operator [](int idx) const {
		return idx < 3 ? v[idx] : 1;
	}
	double& x() {
		return v[0];
	}
	double& y() {
		return v[1];
	}
	double& z() {
		return v[2];
	}
	const double& x() const {
		return v[0];
	}
	const double& y() const {
		return v[1];
	}
	const double& z() const {
		return v[2];
	}
	double lon() const {
		return atan2(y(), x());
	}
	double lat() const {
		return asin(z() / Length());
	}

	point3D operator +(const point3D& t) const {
		point3D ret;
		for (int i = 0; i < 3; ++i) {
			ret.v[i] = v[i] + t.v[i];
		}
		return ret;
	}
	point3D operator -(const point3D& t) const {
		point3D ret;
		for (int i = 0; i < 3; ++i) {
			ret.v[i] = v[i] - t.v[i];
		}
		return ret;
	}
	point3D operator *(const double& t) const {
		point3D ret;
		for (int i = 0; i < 3; ++i) {
			ret.v[i] = v[i] * t;
		}
		return ret;
	}
	point3D operator /(const double& t) const {
		point3D ret;
		for (int i = 0; i < 3; ++i) {
			ret.v[i] = v[i] / t;
		}
		return ret;
	}
	double Length() const {
		double sum = 0;
		for (int i = 0; i < 3; ++i) {
			sum += v[i] * v[i];
		}
		return sqrt(sum);
	}
	double Dot(const point3D& t) const {
		double sum = 0;
		for (int i = 0; i < 3; ++i) {
			sum += v[i] * t.v[i];
		}

		return sum;
	}
	point3D Cross(const point3D& t) const {
		double arr[] = { v[1] * t.v[2] - v[2] * t.v[1], v[2] * t.v[0]
				- v[0] * t.v[2], v[0] * t.v[1] - v[1] * t.v[0] };
		return point3D(arr);
	}
	point3D Normalize() const {
		return point3D(v) / Length();
	}

};

struct sphere {
	point3D c;
	double r;
};

struct plane {
	point3D pnt, norm;
};

struct circle {
	point3D c, norm;
	double r;
};

const int N = 25;
const double earthR =  6370.0;
int n, nn;
double r, lat, lon;
double pntsDist[N * N][N * N];
point3D pnts[N * N];

int dcmp(const double &a, const double &b) {
	if (fabs(a - b) <= EPS)
		return 0;
	return (a > b) * 2 - 1;
}

bool planePlaneIntersection(const plane& p1, const plane& p2, point3D& pnt1,
		point3D& pnt2) {
	point3D dir = p1.norm.Cross(p2.norm);
	if (dcmp(dir.Length(), 0) == 0)
		return 0;
	double n1Dotn1 = p1.norm.Dot(p1.norm);
	double n2Dotn2 = p2.norm.Dot(p2.norm);
	double n1Dotn2 = p1.norm.Dot(p2.norm);
	double d1 = p1.norm.Dot(p1.pnt);
	double d2 = p2.norm.Dot(p2.pnt);

	double determinant = n1Dotn1 * n2Dotn2 - n1Dotn2 * n1Dotn2;
	double c1 = (d1 * n2Dotn2 - d2 * n1Dotn2) / determinant;
	double c2 = (d2 * n1Dotn1 - d1 * n1Dotn2) / determinant;

	pnt1 = p1.norm * c1 + p2.norm * c2;
	pnt2 = pnt1 + dir;
	return 1;
}

bool sphereSphereIntersection(const sphere &s1, const sphere &s2, circle &ret) {
	point3D c1c2 = s2.c - s1.c;
	double d = c1c2.Length();
	if (dcmp(d, s1.r + s2.r) > 0 || dcmp(d, fabs(s1.r - s2.r)) < 0)
		return 0;
	double h = 0.5 + (s1.r * s1.r - s2.r * s2.r) / (2 * d * d);
	ret.c = s1.c + c1c2 * h;
	ret.norm = c1c2.Normalize();
	ret.r = sqrt(s1.r * s1.r - h * h * d * d);
	return 1;
}

vector<point3D> circleCircleIntersectionSamePlane3D(const circle &cr1,
		const circle &cr2) {
	vector<point3D> ret;
	point3D c1c2 = cr2.c - cr1.c;
	double d = c1c2.Length();
	if (dcmp(d, cr1.r + cr2.r) > 0 || dcmp(d, fabs(cr1.r - cr2.r)) < 0)
		return ret;
	double h = 0.5 + (cr1.r * cr1.r - cr2.r * cr2.r) / (2 * d * d);
	double r_i = sqrt(cr1.r * cr1.r - h * h * d * d);
	point3D c_i = cr1.c + c1c2 * h;
	point3D t = c1c2.Cross(cr1.norm).Normalize();
	point3D p_0 = c_i - t * r_i;
	point3D p_1 = c_i + t * r_i;
	ret.push_back(p_0);
	if (dcmp(d, cr1.r + cr2.r) == 0 || dcmp(d, fabs(cr1.r - cr2.r)) == 0)
		return ret;
	ret.push_back(p_1);
	return ret;
}

vector<point3D> sphereCircleIntersection(const sphere &s, const circle &cr) {
	double d = cr.norm.Dot(cr.c - s.c);
	if (dcmp(fabs(d), s.r) > 0)
		return vector<point3D>();
	circle tmp;
	tmp.norm = cr.norm;
	tmp.c = s.c + cr.norm * d;
	tmp.r = sqrt(s.r * s.r - d * d);
	return circleCircleIntersectionSamePlane3D(cr, tmp);
}

double sphereicalDist(const point3D &a, const point3D &b) {
	//return dist(a.lat(), a.lon(), b.lat(), b.lon(), earthR);
	double aa = a.Dot(b) / a.Length() / b.Length();
	aa = max((double) -1.0, min(aa, (double) 1.0));
	aa = acos(aa);
	// 2 * M_PI * erthR: result
	// 2 * M_PI: aa
	return earthR * aa;
}

vector<point3D> intersectLineSphere(const sphere &s, const point3D &p0,
		const point3D &p1) {
	vector<point3D> ret;
	// (p - c).(p - c) = r * r
	// p = p0 + t * (p1 - p0)
	double a, b, c, t1, t2;
	a = (p1 - p0).Dot(p1 - p0);
	b = 2 * (p1 - p0).Dot(p0 - s.c);
	c = (p0 - s.c).Dot(p0 - s.c) - s.r * s.r;
	double det = b * b - 4 * a * c;
	int res;
	if (dcmp(det, 0) < 0)
		return ret;
	else if (dcmp(det, 0.0) == 0)
		det = 0, res = 1;
	else
		res = 2;
	det = sqrt(det);
	t1 = (-b + det) / (2 * a);
	t2 = (-b - det) / (2 * a);
	ret.push_back(p0 + (p1 - p0) * t1);
	if (res == 2)
		ret.push_back(p0 + (p1 - p0) * t2);
	return ret;
}

void addPnts() {
	nn = n;
	//  double v1[] = { 6058.60, -1111.24, 1623.40 };
	//  double v2[] = { 6058.60, 1111.24, 1623.40 };
	//  double v3[] = { 6058.60, 1623.40, -1111.24 };
	//  double v4[] = { 6058.60, 1623.40, 1111.24 };
	//  pnts[nn++] = point3D(v1);
	//  pnts[nn++] = point3D(v2);
	//  pnts[nn++] = point3D(v3);
	//  pnts[nn++] = point3D(v4);
	//  printf("%lf %lf\n", pnts[3].Length(), (pnts[3] - pnts[0]).Length());
	//  printf("%lf %lf\n", pnts[4].Length(), (pnts[4] - pnts[0]).Length());
	//  printf("%lf %lf\n", pnts[5].Length(), (pnts[5] - pnts[0]).Length());
	//  printf("%lf %lf\n", pnts[6].Length(), (pnts[6] - pnts[0]).Length());
	//  return;
	for (int i = 0; i < n; ++i)
		for (int j = i + 1; j < n; ++j) {
			// theta : 2 * M_PI
			// r : 2 * M_PI * erthR
			double theta = r / earthR;
			double h = earthR * cos(theta);
			point3D n1 = pnts[i].Normalize(), n2 = pnts[j].Normalize();
			plane pl1 = { n1 * h, n1 };
			plane pl2 = { n2 * h, n2 };
			point3D ret1, ret2;
			if (!planePlaneIntersection(pl1, pl2, ret1, ret2))
				continue;
			point3D o;
			sphere erth = { o, earthR };

			vector<point3D> retPnts = intersectLineSphere(erth, ret1, ret2);
//			if ((nn == 5 || nn == 9) && retPnts.size())
//				swap(retPnts[0], retPnts[1]);
			for (point3D &retPnt : retPnts) {
				pnts[nn++] = retPnt;
				//        printf("%.2lf %.2lf %.2lf\n", retPnt.x(), retPnt.y(), retPnt.z());
			}
			//      sphere s1 = { pnts[i], r };
			//      sphere s2 = { pnts[j], r };
			//      circle ret;
			//      if (!sphereSphereIntersection(s1, s2, ret))
			//        continue;
			//      sphere erth = { point3D(), erthR };
			//      vector<point3D> retPnts = sphereCircleIntersection(erth, ret);
			//      for (point3D &retPnt : retPnts) {
			//        pnts[nn++] = retPnt;
			//        printf("%.2lf %.2lf %.2lf\n", retPnt.x(), retPnt.y(), retPnt.z());
			//      }
		}
	//  printf("\n");
}

point3D projPointOnLine(const point3D &p0, const point3D &p1,
		const point3D &p) {
	point3D p0p1 = (p1 - p0).Normalize();
	point3D p0p = (p - p0);
	double h = p0p.Dot(p0p1);
	return p0 + p0p1 * h;
}

double computeT(const double &p0, const double &p1,
		const double &p) {
	// p = p0 + T * (p1 - p0)
	return (p - p0) / (p1 - p0);
}

double computeT(const point3D &p0, const point3D &p1, const point3D &p) {
	if (dcmp(p0.x(), p1.x()))
		return computeT(p0.x(), p1.x(), p.x());
	if (dcmp(p0.y(), p1.y()))
		return computeT(p0.y(), p1.y(), p.y());
	return computeT(p0.z(), p1.z(), p.z());
}

struct matrix {
	double arr[4][4];
	matrix operator *(const matrix& m) const {
		matrix ret;
		memset(ret.arr, 0, sizeof(ret.arr));
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				for (int k = 0; k < 4; ++k) {
					ret.arr[i][j] += arr[i][k] * m.arr[k][j];
				}
			}
		}
		return ret;
	}
	point3D operator *(const point3D& m) const {
		point3D ret;
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				ret[i] += arr[i][j] * m[j];
			}
		}
		return ret;
	}
	double& operator()(int i, int j) {
		return arr[i][j];
	}
	const double& operator()(int i, int j) const {
		return arr[i][j];
	}
};
matrix Identity() {
	matrix ret;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			ret(i, j) = i == j;
		}
	}
	return ret;
}
matrix translate(const point3D& v, int dir = 1) {
	matrix ret = Identity();
	for (int i = 0; i < 3; ++i) {
		ret(i, 3) = v[i] * dir;
	}
	return ret;
}
matrix rotateZ(double angle) {
	matrix ret = Identity();
	ret(0, 0) = ret(1, 1) = cos(angle);
	ret(0, 1) = -(ret(1, 0) = sin(angle));
	return ret;
}
matrix transformSystem(const point3D& u, const point3D& v, const point3D& w) {
	matrix ret = Identity();
	for (int j = 0; j < 3; ++j) {
		ret(0, j) = u[j];
		ret(1, j) = v[j];
		ret(2, j) = w[j];
	}
	return ret;
}
matrix ItransformSystem(const point3D& u, const point3D& v, const point3D& w) {
	matrix ret = Identity();
	for (int j = 0; j < 3; ++j) {
		ret(j, 0) = u[j];
		ret(j, 1) = v[j];
		ret(j, 2) = w[j];
	}
	return ret;
}

void getPrep(point3D & w, point3D & v, point3D & u) {
	w = w.Normalize();
	for (int i = 0; i < 3; ++i) {
		if (fabs(w[i]) > EPS) {
			int j = (i + 1) % 3;
			int k = (i + 2) % 3;
			v[i] = w[j];
			v[j] = -w[i];
			v[k] = 0;
			v = v.Normalize();
			break;
		}
	}
	u = v.Cross(w);
}

int circleLineIntersection(const point& p0, const point& p1, const point& cen,
		double rad, point& r1, point & r2) {
	if (same(p0, p1)) {
		if (fabs(lengthSqr(vec(p0, cen)) - (rad * rad)) < EPS) {
			r1 = r2 = p0;
			return 1;
		}
		return 0;
	}
	double a, b, c, t1, t2;
	a = dot(p1 - p0, p1 - p0);
	b = 2 * dot(p1 - p0, p0 - cen);
	c = dot(p0 - cen, p0 - cen) - rad * rad;
	double det = b * b - 4 * a * c;
	int res;
	if (fabs(det) < EPS)
		det = 0, res = 1;
	else if (det < 0)
		res = 0;
	else
		res = 2;
	det = sqrt(det);
	t1 = (-b + det) / (2 * a);
	t2 = (-b - det) / (2 * a);
	r1 = p0 + t1 * (p1 - p0);
	r2 = p0 + t2 * (p1 - p0);
	return res;
}

vector<point3D> getCurveIntersectionWithCircle(point3D p0, point3D p1,
		const point3D &airport) {
	double theta = r / earthR;
	double h = earthR * cos(theta);
	double smallR = earthR * sin(theta);
	point3D w = airport.Normalize(), u, v, circleCenter = w * h;
	getPrep(w, v, u);
	matrix m = transformSystem(u, v, w), im = ItransformSystem(u, v, w);
	p0 = m * p0;
	p1 = m * p1;
	circleCenter = m * circleCenter;
	point P0(p0.x(), p0.y()), P1(p1.x(), p1.y()), ccenter(circleCenter.x(),
			circleCenter.y()), r1, r2;
	int c = circleLineIntersection(P0, P1, ccenter, smallR, r1, r2);
	vector<point3D> ret;
	if (!c)
		return ret;
	double in1[] = { r1.X, r1.Y, circleCenter.z() };
	double in2[] = { r2.X, r2.Y, circleCenter.z() };
	ret.push_back(point3D(in1));
	if (c == 2)
		ret.push_back(point3D(in2));
	for (auto &r : ret)
		r = im * r;
	return ret;
}
map<point3D, int> mp;
double getAngle(const circle &c, point3D p) {
	point3D w = c.norm, u, v, circleCenter = c.c;
	getPrep(w, v, u);
	matrix m = transformSystem(u, v, w);
	p = m * p;
	circleCenter = m * circleCenter;
	double a =  atan2(p.y() , p.x()  );

	return fmod(fmod(a, 2 * M_PI) + 2 * M_PI, 2 * M_PI);
}
bool ok(vector<pair<double, double> >& rngs, const double&st,
		const double&en) {
	double cur = st;
	for (auto rng : rngs) {
		if (dcmp(cur, min(rng.first, en)) < 0)
			return 0;
		cur = max(cur, rng.second);
	}
	return dcmp(cur, en) >= 0;
}
bool canRch(int i, int j) {
	const point3D &p0 = pnts[i], &p1 = pnts[j];
	//point3D p0p1 = (p1 - p0).Normalize();
	circle cr;
	cr.norm = p0.Cross(p1).Normalize();
	cr.r = earthR;
	double st = getAngle(cr, p0);
	double en = getAngle(cr, p1);
	if (en < st)
		swap(st, en);
	if (dcmp(en - st, M_PI) > 0)
		swap(st, en);
	vector<pair<double, double> > rngs;
	for (int k = 0; k < n; ++k) {
		double theta = r / earthR;
		double h = earthR * cos(theta);
		point3D n1 = pnts[k].Normalize();
		plane pl1 = { n1 * h, n1 };
		plane pl2 = { cr.c, cr.norm };
		point3D ret1, ret2;
		if (!planePlaneIntersection(pl1, pl2, ret1, ret2))
			continue;
		sphere erth = { cr.c, cr.r };
		vector<point3D> cur = intersectLineSphere(erth, ret1, ret2);
		//vector<point3D> cur = getCurveIntersectionWithCircle(p0, p1, pnts[k]);
		if (int(cur.size()) < 2)
			continue;
		point3D &a = cur[0], &b = cur[1];
		double st = getAngle(cr, a);
		double en = getAngle(cr, b);
		if (en < st)
			swap(st, en);
		if (dcmp(en - st, M_PI) > 0)
			swap(st, en);
		if (en < st) {
			rngs.push_back(make_pair(st, 2 * M_PI));
			rngs.push_back(make_pair(0, en));
		} else
			rngs.push_back(make_pair(st, en));

//		a = projPointOnLine(p0, p1, a);
//		b = projPointOnLine(p0, p1, b);
//		point3D ab = (b - a).Normalize();
//		if ((p0p1 - ab).Length() > EPS)
//			swap(a, b);
//		double t1 = min(max((double)0.0, computeT(p0, p1, a)),(double) 1.0);
//		double t2 = min(max((double)0.0, computeT(p0, p1, b)), (double)1.0);
//		rngs.push_back(make_pair(t1, t2));
	}
	sort(rngs.begin(), rngs.end());

	if (en < st)
		return ok(rngs, st, 2 * M_PI) && ok(rngs, 0, en);
	return ok(rngs, st, en);
}

void build() {
	for (int i = 0; i < nn; ++i)
		for (int j = i; j < nn; ++j) {
			if (i == j)
				pntsDist[i][j] = 0;
			else
				pntsDist[i][j] = pntsDist[j][i] = (
						canRch(i, j) ?
								sphereicalDist(pnts[i], pnts[j]) : INFINITY);
		}
//	for(auto&a:mp) {int i=a.second,j;
//		for(auto&b:mp)
//			j=b.second,printf("%.3lf ", pntsDist[i][j]);
//		printf("\n");
//	}
//	printf("\n");

	for (int k = 0; k < nn; ++k)
		for (int i = 0; i < nn; ++i)
			for (int j = 0; j < nn; ++j)
				pntsDist[i][j] = min(pntsDist[i][j],
						pntsDist[i][k] + pntsDist[k][j]);

//	for(auto&a:mp) {int i=a.second,j;
//			for(auto&b:mp)
//				j=b.second,printf("%.3lf ", pntsDist[i][j]);
//			printf("\n");
//		}
//		printf("\n");
	//  printf("\n");
}

int notVis[N];
double ds[N];
double dijkstra(int src, int trg, double fuel) {
	fill(ds, ds + N, INFINITY);
	ds[src] = 0;
	int nxt = src;
	for (int i = 0; i < n; i++)
		notVis[i] = i;
	int notVisSz = n;
	do {
		int cur = notVis[nxt];
		notVis[nxt] = notVis[--notVisSz];
		nxt = -1;
		double best = INFINITY;
		for (int k = 0; k < notVisSz; k++) {
			int j = notVis[k];
			double dj = ds[cur] + pntsDist[cur][j];
			if (dcmp(pntsDist[cur][j], fuel) <= 0 && ds[j] > dj)
				ds[j] = dj;
			if (best > ds[j])
				best = ds[j], nxt = k;
		}
	} while (nxt != -1);
	return ds[trg];
}

//void test() {
//  double v1[] = { 3, 4, 2 };
//  double v2[] = { 5, 2, 6 };
//  double v3[] = { 4, 4, 2 };
//  double v4[] = { 3, 5, 2 };
//  double v5[] = { 3, 4, 3 };
//  double v6[] = { 4, 3, 3 };
//
//  point3D p1(v1), p2(v2), p3(v3), p4(v4), p5(v5), p6(v6);
//  point3D n1 = (p3 - p1).Cross(p2 - p1);
//  point3D n2 = (p6 - p4).Cross(p5 - p4);
//  plane P1 = { p1, n1 };
//  plane P2 = { p4, n2 };
//
//  point3D ret1, ret2;
//  if (!planePlaneIntersection(P1, P2, ret1, ret2))
//    return;
//  printf("%.2lf %.2lf %.2lf\n", ret1.x(), ret1.y(), ret1.z());
//  printf("%.2lf %.2lf %.2lf\n", ret2.x(), ret2.y(), ret2.z());
//  double t = computeT(ret1, ret2, p3);
//  point3D tmp = ret1 + (ret2 - ret1) * t;
//  printf("%.2lf %.2lf %.2lf\n", tmp.x(), tmp.y(), tmp.z());
//}

bool operator <(const point3D &a, const point3D &b) {
	for (int i = 0; i < 3; ++i)
		if (fabs(a[i] - b[i]) > EPS)
			return a[i] < b[i];
	return false;
}

int main() {
	// freopen("test.in", "rt", stdin);
	//freopen("wa.txt", "wt", stdout);
	int tst = 1;
	double tmp, tmp1;
	while (~scanf("%d%lf", &n, &tmp)) {
		r = tmp;
		printf("Case %d:\n", tst++);
		for (int i = 0; i < n; ++i) {
			scanf("%lf%lf", &tmp, &tmp1);
			lon = tmp, lat = tmp1;
			pnts[i] = point3D((lon * M_PI) / 180.0, (lat * M_PI) / 180.0,
					earthR);
		}
		addPnts();
		for (int i = 0; i < nn; i++) {
			mp[pnts[i]] = i;
		}
		build();
//		for (int i = 0; i < nn; i++) {
//			printf("%d %.3lf %.3lf %.3lf\n", i, pnts[i][0], pnts[i][1] ,
//					pnts[i][2]);
//		}
		int q, src, trg;
		double fuel;
		scanf("%d", &q);

		while (q--) {
			scanf("%d%d%lf", &src, &trg, &tmp);
			fuel = tmp;

			double ret = dijkstra(src - 1, trg - 1, fuel);
			if (ret == INFINITY)
				printf("impossible\n");
			else
				printf("%.3lf\n", (double) (ret ));
		}
	}
	return 0;
}
