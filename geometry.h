#include <iostream>
#include <cstring>
#include <vector>
#include <cmath>

struct Point {
  double x = 0;
  double y = 0;

  Point(double x, double y) : x(x), y(y) {}

  Point() = default;
};

struct Vector {
  Point point = Point();

  Vector(const Point p) : point(p) {}

  Vector(const Point p1, const Point p2) : point(Point(p2.x - p1.x, p2.y - p1.y)) {}

  double squaredLength() const {
    return (point.x * point.x + point.y * point.y);
  }

  double length() const {
    return sqrt(squaredLength());
  }

  Point getPoint() {
    return point;
  }

  Vector &operator=(Vector v) {
    point = v.point;
    return *this;
  }

  Vector operator-() const {
    Point p = point;
    p.x = -p.x;
    p.y = -p.y;
    return p;
  }

  Vector RotateLeft() {
    Point temp;
    temp.x = -point.y;
    temp.y = point.x;
    return temp;
  }

  void RotateLeft(double angle) {
    double rad = angle / (180 * M_PI);
    double t = point.x;
    point.x = point.x * cos(rad) - point.y * sin(rad);
    point.y = t * sin(rad) + point.y * cos(rad);
  }

};


Vector operator+(Vector v1, Vector v2) {
  v1.point.x = (v1.point.x + v2.point.x);
  v1.point.y = (v1.point.y + v2.point.y);
  return v1;
}

Vector operator-(Vector v1, Vector v2) {
  v2 = -v2;
  return v1 + v2;
}

double dist(Point p1, Point p2) {
  auto v1 = Vector(p1);
  auto v2 = Vector(p2);
  return (v1 - v2).length();
}

Vector operator/(Vector v1, double div) {
  v1.point.x /= div;
  v1.point.y /= div;
  return v1;
}

Vector operator*(Vector v, double mul) {
  v.point.x *= mul;
  v.point.y *= mul;
  return v;
}

Vector operator*(double mul, Vector v) {
  v.point.x *= mul;
  v.point.y *= mul;
  return v;
}

Vector operator/(double div, Vector v) {
  v.point.x /= div;
  v.point.y /= div;
  return v;
}

Point mid(Point p1, Point p2) {
  auto v1 = Vector(p1);
  auto v2 = Vector(p2);
  return ((v1 + v2) / 2).getPoint();
}

double koef(Vector v1, Vector v2) {
  return v1.squaredLength() / v2.squaredLength();
}

double scalar(Vector v1, Vector v2) {
  return v1.point.x * v2.point.x + v1.point.y * v2.point.y;
}

Point RotatePoint(Point p, Point center, double angle) {
  auto v = Vector(center, p);
  v.RotateLeft(angle);
  v = v + center;
  return v.getPoint();
}

double angle(Vector v1, Vector v2) {
  return acos((scalar(v1, v2)) / (v1.length() * v2.length()));
}

bool operator==(const Point p1, const Point p2) {
  return (fabs(p1.x - p2.x) < 1e-4 && fabs(p1.y - p1.y) < 1e-4);
}

bool operator!=(const Point p1, const Point p2) {
  return !(p1 == p2);
}

bool operator==(const Vector v1, const Vector v2) {
  return v1.point == v2.point;
}

bool operator!=(const Vector v1, const Vector v2) {
  return !(v1 == v2);
}

class Line {
 public:
  double a_ = 1;
  double b_ = 0;
  double c_ = 0;

  Line() = default;

  Line(double a, double b, double c) : a_(a), b_(b), c_(c) {}

  Line(double k, double b) : a_(-k), b_(1), c_(-b) {}

  Line(double k, Point p) : a_(k), b_(1), c_(-k * p.x - p.y) {}

  Line(Point p, double k) : a_(k), b_(1), c_(-k * p.x - p.y) {}

  Line(Point p1, Point p2) {
    if (fabs(p1.x - p2.x) < 1e-4) {
      a_ = 1;
      b_ = 0;
      c_ = -p1.x;
    } else {
      a_ = (p1.y - p2.y) / (p2.x - p1.x);
      b_ = 1;
      c_ = -(a_ * p1.x + p1.y);
    }
  }

  Line getPerpend(Point p) const {
    double ta = -b_;
    double tb = a_;
    double tc = p.x * b_ - p.y * a_;
    if (abs(tb) > 1e-4) {
      ta /= tb;
      tc /= tb;
      tb = 1;
    } else {
      tb = 0;
    }
    Line l(ta, tb, tc);
    return l;

  }

};

Point intersection(Line l1, Line l2) {
  double det = l1.a_ * l2.b_ - l2.a_ * l1.b_;
  Point p;
  p.x = (l1.b_ * l2.c_ - l2.b_ * l1.c_) / det;
  p.y = (l2.a_ * l1.c_ - l1.a_ * l2.c_) / det;
  return p;
}

Point ReflexPoint(Point p, Line l) {
  Line l1 = l.getPerpend(p);
  Point inter = intersection(l1, l);
  auto v = Vector(p, inter);
  v = 2 * v;
  v = v + Vector(p);
  return v.getPoint();
}

Point ScalePoint(Point p, Point center, double coef) {
  Vector v(center, p);
  v = center + (v * coef);
  return v.getPoint();
}


bool operator==(const Line l1, const Line l2) {
  if (l1.b_ == l2.b_ && fabs(l1.a_ - l2.a_) < 1e-4) {
    return true;
  }
  return false;
}

bool operator!=(const Line l1, const Line l2) {
  return !(l1 == l2);
}

class Shape {
 public:
  virtual double perimeter() const = 0;

  virtual double area() const = 0;

  virtual bool isCongruentTo(const Shape &another) const = 0;

  virtual bool isSimilarTo(const Shape &another) const = 0;

  virtual bool containsPoint(Point p) const = 0;

  virtual bool operator==(const Shape &another) const = 0;

  virtual bool operator!=(const Shape &another) const = 0;

  virtual void scale(Point, double) = 0;

  virtual void reflex(Line) = 0;

  virtual void reflex(Point) = 0;

  virtual void rotate(Point, double) = 0;

  virtual ~Shape() = 0;
};

Shape::~Shape() {}

class Polygon : public Shape {
 public:
  Polygon() = default;

  Polygon(const std::vector<Point> points) : points_(points) {}

  size_t verticesCount() const {
    return points_.size();
  }

  std::vector<Point> getVertices() const {
    return points_;
  }

  double perimeter() const override {
    double ans = 0;
    size_t sz = verticesCount();
    for (size_t i = 0; i < sz; ++i) {
      ans += (Vector(points_[i], points_[(i + 1) % sz]).length());
    }
    return ans;

  }

  double area() const override {
    double ans = 0;
    size_t sz = verticesCount();
    for (size_t i = 0; i < sz; ++i) {
      ans += (points_[i].x + points_[(i + 1) % sz].x) * (points_[i].y - points_[(i + 1) % sz].y);
    }
    return fabs(ans / 2);
  }

  bool isConvex() const {
    size_t sz = verticesCount();
    size_t cnt = 0;
    for (size_t i = 0; i < sz; ++i) {
      Vector v1 = Vector(points_[i], points_[(sz + i - 1) % sz]);
      Vector v2 = Vector(points_[i], points_[(i + 1) % sz]);
      if (v1.point.x * v2.point.y > v1.point.y * v2.point.x) {
        ++cnt;
      }
    }
    if (cnt == 0 || cnt == sz) return true;
    return false;
  }

  bool operator==(const Shape &another) const override {
    auto *second = dynamic_cast<const Polygon *>(&another);
    if (second == nullptr) {
      return false;
    }
    if (verticesCount() != second->verticesCount()) {
      return false;
    }
    size_t sz = verticesCount();
    std::vector<Point> v1 = getVertices();
    std::vector<Point> v2 = second->getVertices();
    for (size_t i = 0; i < sz; ++i) {
      bool flag = true;
      for (size_t j = 0; j < sz; ++j) {
        if (v1[(i + j) % sz] != v2[j]) {
          flag = false;
          break;
        }
      }
      if (flag) {
        return true;
      }
      flag = true;
      for (size_t j = 0; j < sz; ++j) {
        if (v1[(sz + i - j) % sz] != v2[j]) {
          flag = false;
          break;
        }
      }
      if (flag) {
        return true;
      }
    }
    return false;
  }

  bool operator!=(const Shape &another) const override {
    return !(*this == another);
  }

  bool isSimilarTo(const Shape &another) const override {
    auto *second = dynamic_cast<const Polygon *>(&another);
    if (second == nullptr) {
      return false;
    }
    if (verticesCount() != second->verticesCount()) {
      return false;
    }
    size_t sz = verticesCount();
    std::vector<Vector> v1;
    for (size_t i = 1; i < sz; ++i) {
      v1.push_back(Vector(points_[0], points_[i]));
    }
    for (size_t i = 0; i < sz; ++i) {
      std::vector<Vector> v2;
      v2.clear();
      for (size_t j = 1; j < sz; ++j) {
        v2.push_back(Vector(points_[i], points_[(i + j) % sz]));
      }
      double k = koef(v1[0], v2[0]);
      bool flag = true;
      for (size_t j = 1; j < sz; ++j) {
        if (fabs(k - koef(v1[j], v2[j])) > 1e-4) {
          flag = false;
          break;
        }
      }
      if (flag) return true;
      k = koef(v1[0], v2[sz - 1]);
      flag = true;
      for (size_t j = 1; j < sz; ++j) {
        if (fabs(k - koef(v1[j], v2[sz - j - 1])) > 1e-4) {
          flag = false;
          break;
        }
      }
      if (flag) return true;
    }
    return false;
  }

  bool isCongruentTo(const Shape &another) const override {
    if (!isSimilarTo(another)) return false;
    auto *second = dynamic_cast<const Polygon *>(&another);
    if (fabs(area() - second->area()) < 1e-4) {
      return true;
    }
    return false;
  }

  bool containsPoint(Point point) const override {
    double sum = 0;
    size_t sz = verticesCount();
    for (size_t i = 0; i < sz; ++i) {
      auto v1 = Vector(point, points_[i]);
      auto v2 = Vector(point, points_[(i + 1) % sz]);
      sum += angle(v1, v2);
    }
    if (fabs(sum - 2 * M_PI) < 1e-4) return true;
    return false;
  }

  void rotate(Point center, double angle) override {
    for (size_t i = 0; i < verticesCount(); ++i) {
      points_[i] = RotatePoint(points_[i], center, angle);
    }
  }

  void reflex(Point center) override {
    rotate(center, 180);
  }

  void reflex(Line l) override {
    for (size_t i = 0; i < verticesCount(); ++i) {
      points_[i] = ReflexPoint(points_[i], l);
    }
  }

  void scale(Point center, double coefficient) override {
    for (size_t i = 0; i < verticesCount(); ++i) {
      points_[i] = ScalePoint(points_[i], center, coefficient);
    }
  }

 protected:
  std::vector<Point> points_ = std::vector<Point>();
};

class Ellipse : public Shape {
 public:
  Ellipse(Point p1, Point p2, double a) : first_focus_(p1), second_focus_(p2), a_(a / 2) {}

  Ellipse() = default;

  std::pair<Point, Point> focuses() const {
    return std::make_pair(first_focus_, second_focus_);
  }

  double eccentricity() const {
    auto v = Vector(first_focus_, center());
    return v.length() / (a_);
  }

  double getb() const {
    double c = dist(first_focus_, center());
    return sqrt(a_ * a_ - c * c);
  }

  virtual Point center() const {
    return mid(first_focus_, second_focus_);
  }

  std::pair<Line, Line> directrices() const {
    auto v_center = Vector(center());
    auto v_focus = Vector(first_focus_, center());
    v_focus = v_focus * (a_ * a_) / v_focus.squaredLength();
    auto b = v_focus.RotateLeft();
    auto direct = v_center + v_focus;
    Line first_line = Line((direct).getPoint(), (direct + b).getPoint());
    direct = v_center - v_focus;
    Line second_line = Line((direct).getPoint(), (direct + b).getPoint());
    return std::make_pair(first_line, second_line);
  }

  double area() const override {
    double b = getb();
    return M_PI * b * a_;
  }

  double perimeter() const override {
    double b = getb();
    double ans = M_PI * (3 * (a_ + b) - sqrt((3 * a_ + b) * (3 * b + a_)));
    return ans;
  }

  bool operator!=(const Shape &another) const override {
    return !(*this == another);
  }

  bool operator==(const Shape &another) const override {
    auto *second = dynamic_cast<const Ellipse *>(&another);
    if (second == nullptr) {
      return false;
    }
    if (fabs(a_ - second->a_) < 1e-4) {
      if (first_focus_ == second->first_focus_ && second_focus_ == second->second_focus_) {
        return true;
      }
      if (first_focus_ == second->second_focus_ && second_focus_ == second->first_focus_) {
        return true;
      }
    }
    return false;
  }

  bool isCongruentTo(const Shape &another) const override {
    auto *second = dynamic_cast<const Ellipse *>(&another);
    if (second == nullptr) {
      return false;
    }
    if (fabs(eccentricity() - second->eccentricity()) > 1e-4) {
      return false;
    }
    if (fabs(area() - second->area()) < 1e-4) {
      return true;
    }
    return false;
  }

  bool isSimilarTo(const Shape &another) const override {
    auto *second = dynamic_cast<const Ellipse *>(&another);
    if (second == nullptr) {
      return false;
    }
    if (fabs(eccentricity() - second->eccentricity()) < 1e-4) {
      return true;
    }
    return false;
  }

  bool containsPoint(Point point) const override {
    std::pair<Point, Point> pr = focuses();
    double sum = dist(point, pr.first) + dist(point, pr.second);
    if (sum < 2 * a_) {
      return true;
    }
    return false;
  }

  void rotate(Point center, double angle) override {
    first_focus_ = RotatePoint(first_focus_, center, angle);
    second_focus_ = RotatePoint(second_focus_, center, angle);
  }

  void reflex(Point center) override {
    rotate(center, 180);
  }

  void reflex(Line l) override {
    first_focus_ = ReflexPoint(first_focus_, l);
    second_focus_ = ReflexPoint(second_focus_, l);
  }

  void scale(Point center, double coefficient) override {
    first_focus_ = ScalePoint(first_focus_, center, coefficient);
    second_focus_ = ScalePoint(second_focus_, center, coefficient);
    a_ = coefficient * a_;
  }

 protected:
  Point first_focus_;
  Point second_focus_;
  double a_{0};
};

class Circle : public Ellipse {
 public:
  Circle(Point p, double r) {
    first_focus_ = p;
    second_focus_ = p;
    a_ = r;
  }

  double radius() const {
    return a_;
  }

  Point center() const override {
    return first_focus_;
  }

  double area() const override {
    return M_PI * radius() * radius();
  }

  double perimeter() const override {
    double ans = 2 * M_PI * radius();
    return ans;
  }


};

class Rectangle : public Polygon {
 public:
  Rectangle(Point p1, Point p3, double coefficient) {
    if (coefficient < 1) {
      coefficient = 1 / coefficient;
    }
    auto diameter = Vector(p1, p3);
    double b_in_sqr = diameter.squaredLength() / (1 + coefficient * coefficient);
    coefficient = 1 / coefficient;
    double a_in_sqr = diameter.squaredLength() / (1 + coefficient * coefficient);
    double h = sqrt(a_in_sqr * b_in_sqr / diameter.squaredLength());
    Vector d_b = diameter * (b_in_sqr / diameter.squaredLength());
    Vector vector_h = (diameter.RotateLeft()) * (h / diameter.length());
    Point p2 = (Vector(p1) + d_b + vector_h).getPoint();
    Line l1 = Line(p2, mid(p1, p3));
    Line l2 = Line(p1, p2);
    Line l3 = l2.getPerpend(p1);
    Point p4 = intersection(l1, l3);
    points_.clear();
    points_.push_back(p1);
    points_.push_back(p2);
    points_.push_back(p3);
    points_.push_back(p4);
  }

  Point center() const {
    return ((Vector(points_[0]) + Vector(points_[3])) / 2).getPoint();
  }
};

class Square : public Rectangle {
 public:
  Square(Point p1, Point p2) : Rectangle(p1, p2, 1) {}

  Circle circumscribedCircle() {
    Point cir = center();
    double r = dist(points_[0], cir);
    return Circle(cir, r);
  }

  Circle inscribedCircle() {
    Point cir_center = center();
    Point ma = mid(points_[0], points_[1]);
    double r = dist(ma, cir_center);
    return Circle(cir_center, r);
  }
};

class Triangle : public Polygon {
 public:
  Triangle(Point p1, Point p2, Point p3) {
    points_.clear();
    points_.push_back(p1);
    points_.push_back(p2);
    points_.push_back(p3);
  }

  Circle ninePointsCircle() const {
    Point middle_a = mid(points_[1], points_[2]);
    Point middle_b = mid(points_[0], points_[2]);
    Point middle_c = mid(points_[1], points_[0]);
    Triangle med_tr = Triangle(middle_a, middle_b, middle_c);
    return med_tr.circumscribedCircle();
  }

  Line EulerLine() const {
    Point h = orthocenter();
    Point m = centroid();
    return Line(h, m);
  }

  Point orthocenter() const {
    Line l1(points_[0], points_[1]);
    Line l2(points_[0], points_[2]);
    Line h1 = l1.getPerpend(points_[2]);
    Line h2 = l2.getPerpend((points_[1]));
    return intersection(h1, h2);
  }

  Point centroid() const {
    auto v1 = Vector(points_[0]);
    auto v2 = Vector(points_[1]);
    auto v3 = Vector(points_[2]);
    return ((v1 + v2 + v3) / 3).getPoint();
  }

  Circle inscribedCircle() const {
    auto a = Vector(points_[2], points_[1]);
    auto b = Vector(points_[2], points_[0]);
    auto c = Vector(points_[0], points_[1]);
    auto C = Vector(points_[2]);
    double v_len = perimeter() / 2 - c.length();
    a = a * (v_len / a.length());
    b = b * (v_len / b.length());
    auto la = Line(points_[1], points_[2]);
    auto lb = Line(points_[2], points_[0]);
    auto l1 = la.getPerpend((C + a).getPoint());
    auto l2 = lb.getPerpend((C + b).getPoint());
    Point center = intersection(l1, l2);
    double r = (C + a - Vector(center)).length();
    return Circle(center, r);
  }

  Circle circumscribedCircle() const {
    auto a = Line(points_[1], points_[2]);
    auto b = Line(points_[2], points_[0]);
    auto bm = mid(points_[0], points_[2]);
    auto am = mid(points_[1], points_[2]);
    Line l1 = a.getPerpend(am);
    Line l2 = b.getPerpend(bm);
    Point center = intersection(l1, l2);
    double radius = dist(center, points_[0]);
    return Circle(center, radius);
  }
};



