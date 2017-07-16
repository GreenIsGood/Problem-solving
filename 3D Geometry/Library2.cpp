
#define EPS 1e-9

double ONE = 1;
struct point3D {
        double v[3];
        point3D() {
                for (int i = 0; i < 3; ++i) {
                        this->v[i] = 0;
                }
        }
        point3D(double v[3]) {
                for (int i = 0; i < 3; ++i) {
                        this->v[i] = v[i];
                }
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
        double Length() {
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
                double arr[] = { v[1] * t.v[2] - v[2] * t.v[1], v[2] * t.v[0] - v[0]
                                * t.v[2], v[0] * t.v[1] - v[1] * t.v[0] };
                return point3D(arr);
        }
        point3D Normalize() {
                return point3D(v) / Length();
        }
};
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
matrix rotate(const point3D& p, const point3D& q, double angle) {
        point3D w((q - p).Normalize()), u, v;
        getPrep(w, v, u);

        return translate(p, 1) * ItransformSystem(u, v, w) * rotateZ(angle)
                        * transformSystem(u, v, w) * translate(p, -1);
}
bool linePlaneIntersect(const point3D& p, const point3D& q, const point3D& pp,
                const point3D& N, point3D& ret) {
        double d = (q - p).Dot(N);
        if (fabs(d) < EPS)
                return false;

        double t = (pp - p).Dot(N) / d;
        ret = p + (q - p) * t;
        return true;
}
point3D tetra_center(const point3D & a, const point3D & b, const point3D & c,
                const point3D & d) {
        return (a + b + c + d) / 4;
}
double tetra_volume(const point3D & a, const point3D & b, const point3D & c,
                const point3D & d) {
        return fabs((a - d).Dot((b - d).Cross(c - d))) / 6;
}
