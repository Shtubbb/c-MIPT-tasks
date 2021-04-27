#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <math.h>

class BigInteger {
 public:
  static BigInteger BinPow(BigInteger num, size_t degree) {
    BigInteger res = 1;
    while (degree) {
      if (degree % 2 == 1) {
        res *= num;
      }
      num *= num;
      degree /= 2;
    }
    return res;
  }

  void MakePos() {
    sign_ = true;
  }


  BigInteger() = default;

  BigInteger(int val) {
    sign_ = true;
    if (val < 0) {
      sign_ = false;
      val = -val;
    }
    if (val < base_) {
      num_[0] = val;
    } else {
      num_[0] = (val % base_);
      num_.push_back(val / base_);
    }
  }

  BigInteger(std::string str) {
    sign_ = true;
    int length = str.length();
    if (str[0] == '+') {
      str = str.substr(1, length - 1);
    }
    if (str[0] == '-') {
      sign_ = false;
      str = str.substr(1, length - 1);
    }
    if (str[0] == '0') {
      *this = 0;
      return;
    }
    length = str.length();
    int temp;
    if (!num_.empty()) {
      num_.clear();
    }
    for (int i = length; i >= 0; i -= 9) {
      temp = 0;
      if (i < 9) {
        for (int j = 0; j < i; ++j) {
          temp = temp * 10 + str[j] - '0';
        }
      } else {
        for (int j = i - 9; j < i; ++j) {
          temp = temp * 10 + str[j] - '0';
        }
      }
      num_.push_back(temp);
    }
  }

  BigInteger(const BigInteger &val) = default;

  bool GetSign() const {
    return sign_;
  }

  std::vector<int> GetNum() const {
    return num_;
  }

  int operator[](size_t i) const {
    return num_[i];
  }

  size_t GetNumSize() {
    return num_.size();
  }

  BigInteger &operator=(BigInteger val) {
    std::swap(sign_, val.sign_);
    std::swap(num_, val.num_);
    return *this;
  }

  BigInteger operator-() const {
    BigInteger result = *this;
    result.sign_ = !result.sign_;
    result.Fix();
    return result;
  }

  explicit operator bool() const {
    return !(num_.size() == 1 && num_[0] == 0);
  }

  explicit operator double() const {
    double ans = 0;
    double base = 1;
    size_t size = num_.size();
    for (size_t i = 0; i < size; ++i) {
      ans += static_cast<double>(num_[i]) * base;
      base *= base_;
    }
    if (!sign_) {
      ans = -ans;
    }
    return ans;
  }

  BigInteger &operator+=(const BigInteger &val) {
    if (!sign_ && !val.sign_) {
      BigInteger temp = -*this;
      temp += -val;
      return *this = -temp;
    }
    if (!sign_) {
      BigInteger temp = val;
      temp -= -*this;
      return *this = temp;
    }
    if (!val.sign_) {
      return (*this -= -val);
    }
    for (size_t i = 0; i < val.num_.size(); ++i) {
      if (i >= num_.size()) {
        num_.push_back(0);
      }
      num_[i] += val.num_[i];
    }
    for (size_t i = 0; i < num_.size(); ++i) {
      if (num_[i] >= base_) {
        num_[i] -= base_;
        if (i == num_.size() - 1) {
          num_.push_back(0);
        }
        num_[i + 1] += 1;
      }
    }
    Fix();
    return *this;
  }

  BigInteger &operator*=(const BigInteger &val) {
    BigInteger temp;
    temp.num_.resize(val.num_.size() + num_.size(), 0);
    for (size_t i = 0; i < val.num_.size(); ++i) {
      if (val.num_[i] == 0) {
        continue;
      }
      int next = 0;
      size_t j;
      for (j = 0; j < num_.size(); ++j) {
        auto temp_num = static_cast<long long>(val.num_[i]) * num_[j] + next + temp.num_[i + j];
        temp.num_[i + j] = temp_num % base_;
        next = temp_num / base_;
      }
      while (next) {
        long long temp_num = static_cast<long long>(temp.num_[i + j]) + next;
        temp.num_[i + j] = temp_num % base_;
        next = temp_num / base_;
        ++j;
      }
    }
    temp.Fix();
    if (val.sign_ == sign_) {
      return *this = temp;
    }
    return *this = -temp;
  }

  BigInteger &operator/=(const BigInteger &);

  BigInteger &operator-=(const BigInteger &);

  BigInteger &operator%=(const BigInteger &);

  BigInteger &operator--() {
    *this -= 1;
    return *this;
  }

  BigInteger &operator++() {
    *this += 1;
    return *this;
  }

  BigInteger operator--(int) {
    BigInteger temp = *this;
    *this -= 1;
    return temp;
  }

  BigInteger operator++(int) {
    BigInteger temp = *this;
    *this += 1;
    return temp;
  }

  std::string toString() const {
    std::string result;
    if (!sign_) {
      result = "-";
    }
    if (!num_.empty()) {
      result += std::to_string(num_.back());
    }
    for (int i = static_cast<int>(num_.size()) - 2; i >= 0; --i) {
      int temp_base = base_;
      int temp_val = num_[i];
      for (size_t j = 0; j < 9; ++j) {
        temp_base /= 10;
        result += '0' + temp_val / temp_base;
        temp_val %= temp_base;
      }
    }
    return result;
  }

 private:
  static const int base_ = 1'000'000'000;
  std::vector<int> num_{0};
  bool sign_{true};

  friend std::istream &operator>>(std::istream &, BigInteger &);

  void Fix() {
    if (num_.empty()) {
      num_.push_back(0);
    }
    while (num_.size() > 1 && num_.back() == 0) {
      num_.pop_back();
    }
    if (num_.size() == 1 && num_[0] == 0) {
      sign_ = true;
    }
  }

};

bool operator==(const BigInteger &first_val, const BigInteger &second_val) {
  return (first_val.GetSign() == second_val.GetSign() && first_val.GetNum() == second_val.GetNum());
}

bool operator!=(const BigInteger &first_val, const BigInteger &second_val) {
  return !(first_val == second_val);
}

bool operator>(const BigInteger &first_val, const BigInteger &second_val) {
  bool first_sign = first_val.GetSign();
  bool second_sign = second_val.GetSign();
  if (first_sign && !second_sign) {
    return true;
  }
  if (!first_sign && second_sign) {
    return false;
  }
  if (!first_sign && !second_sign) {
    return -second_val > -first_val;
  }
  std::vector<int> first_num = first_val.GetNum();
  std::vector<int> second_num = second_val.GetNum();
  if (first_num.size() > second_num.size()) {
    return true;
  }
  if (first_num.size() < second_num.size()) {
    return false;
  }
  for (int i = static_cast<int>(first_num.size()) - 1; i >= 0; --i) {
    if (first_num[i] < second_num[i]) {
      return false;
    }
    if (first_num[i] > second_num[i]) {
      return true;
    }
  }
  return false;

}

bool operator<(const BigInteger &first_val, const BigInteger &second_val) {
  bool first_sign = first_val.GetSign();
  bool second_sign = second_val.GetSign();
  if (!first_sign && second_sign) {
    return true;
  }
  if (first_sign && !second_sign) {
    return false;
  }
  if (!first_sign && !second_sign) {
    return -second_val < -first_val;
  }
  std::vector<int> first_num = first_val.GetNum();
  std::vector<int> second_num = second_val.GetNum();
  if (first_num.size() < second_num.size()) {
    return true;
  }
  if (first_num.size() > second_num.size()) {
    return false;
  }
  for (int i = static_cast<int>(first_num.size()) - 1; i >= 0; --i) {
    if (first_num[i] > second_num[i]) {
      return false;
    }
    if (first_num[i] < second_num[i]) {
      return true;
    }

  }
  return false;
}

bool operator>=(const BigInteger &first_val, const BigInteger &second_val) {
  return !(first_val < second_val);
}

bool operator<=(const BigInteger &first_val, const BigInteger &second_val) {
  return !(first_val > second_val);
}

BigInteger operator+(const BigInteger &first_val, const BigInteger &second_val) {
  BigInteger temp = first_val;
  temp += second_val;
  return temp;
}

BigInteger operator-(const BigInteger &first_val, const BigInteger &second_val) {
  BigInteger temp = first_val;
  temp -= second_val;
  return temp;
}

BigInteger operator*(const BigInteger &first_val, const BigInteger &second_val) {
  BigInteger temp = first_val;
  temp *= second_val;
  return temp;
}

BigInteger operator/(const BigInteger &first_val, const BigInteger &second_val) {
  BigInteger temp = first_val;
  temp /= second_val;
  return temp;
}

BigInteger operator%(const BigInteger &first_val, const BigInteger &second_val) {
  BigInteger temp = first_val;
  temp %= second_val;
  return temp;
}

BigInteger &BigInteger::operator-=(const BigInteger &val) {
  if (!GetSign() && !val.GetSign()) {
    BigInteger temp = -val;
    temp -= -*this;
    return *this = temp;
  }
  if (!GetSign()) {
    BigInteger temp = val;
    temp += -*this;
    return *this = -temp;
  }
  if (!val.GetSign()) {
    return (*this += -val);
  }
  if (*this < val) {
    BigInteger temp = val;
    temp -= *this;
    return *this = -temp;
  }
  for (size_t i = 0; i < val.num_.size(); ++i) {
    num_[i] -= val.num_[i];
  }
  for (size_t i = 0; i < num_.size(); ++i) {
    if (num_[i] < 0) {
      num_[i] += base_;
      num_[i + 1] -= 1;
    }
  }
  Fix();
  return *this;
}

BigInteger &BigInteger::operator/=(const BigInteger &val) {
  BigInteger ans;
  if (sign_ != val.sign_) {
    ans.sign_ = false;
  }
  if (val.num_.size() == 1 && num_.size() == 1) {
    ans.num_[0] = num_[0] / val.num_[0];
    *this = ans;
    return *this;
  }
  BigInteger front = 0;
  ans.num_.clear();
  front.num_.clear();
  for (int i = static_cast<int>(num_.size()) - 1; i >= 0; --i) {
    for (size_t j = 0; j < front.num_.size() / 2; ++j) {
      std::swap(front.num_[j], front.num_[front.num_.size() - j - 1]);
    }
    front.num_.push_back(num_[i]);
    for (size_t j = 0; j < front.num_.size() / 2; ++j) {
      std::swap(front.num_[j], front.num_[front.num_.size() - j - 1]);
    }
    front.Fix();
    BigInteger cur;
    if (front >= val) {
      int l = 1, r = base_;
      while (r - l > 1) {
        int mid = (r + l) / 2;
        cur = val * mid;
        if (cur <= front) {
          l = mid;
        } else {
          r = mid;
        }
      }
      ans.num_.push_back(l);
      cur = val;
      cur *= l;
      front -= cur;
    } else {
      ans.num_.push_back(0);
    }
  }
  for (size_t j = 0; j < ans.num_.size() / 2; ++j) {
    std::swap(ans.num_[j], ans.num_[ans.num_.size() - j - 1]);
  }
  *this = ans;
  Fix();
  return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &val) {
  return *this = (*this - (*this / val) * val);
}

BigInteger operator ""_bi(unsigned long long val) {
  return BigInteger(val);
}

std::ostream &operator<<(std::ostream &out, const BigInteger &num) {
  std::string string_num = num.toString();
  out << string_num;
  return out;
}

std::istream &operator>>(std::istream &in, BigInteger &val) {
  std::string str;
  in >> str;
  val = str;
  return in;
}

BigInteger Gcd(BigInteger first, BigInteger second) {
  first.MakePos();
  second.MakePos();
  while (second) {
    first %= second;
    std::swap(first, second);
  }
  return first;
}


class Rational {
 public:

  friend std::istream &operator>>(std::istream &, Rational &);

  Rational() = default;

  Rational(const BigInteger &val) : n_(val), m_(1) {}

  Rational(const BigInteger &val_first, const BigInteger &val_second) : n_(val_first), m_(val_second) {}

  Rational(int val) : n_(val), m_(1) {}

  Rational(const Rational &) = default;

  Rational operator-() const {
    Rational temp = *this;
    temp.n_ = -temp.n_;
    return temp;
  }

  BigInteger GetN() const {
    return n_;
  }

  BigInteger GetM() const {
    return m_;
  }


  Rational &operator=(Rational val) {
    std::swap(val.n_, n_);
    std::swap(val.m_, m_);
    return *this;
  }

  Rational &operator+=(const Rational &val) {
    n_ = n_ * val.m_ + m_ * val.n_;
    m_ *= val.m_;
    BigInteger gcd = Gcd(n_, m_);
    m_ /= gcd;
    n_ /= gcd;
    FixRational();
    return *this;
  }

  Rational &operator-=(const Rational &val) {
    n_ = n_ * val.m_ - m_ * val.n_;
    m_ *= val.m_;
    BigInteger gcd = Gcd(n_, m_);
    m_ /= gcd;
    n_ /= gcd;
    FixRational();
    return *this;
  }

  Rational &operator*=(Rational val) {
    BigInteger gcd = Gcd(n_, val.m_);
    val.m_ /= gcd;
    n_ /= gcd;
    gcd = Gcd(m_, val.n_);
    val.n_ /= gcd;
    m_ /= gcd;
    n_ *= val.n_;
    m_ *= val.m_;
    FixRational();
    return *this;
  }

  Rational &operator/=(Rational val) {
    BigInteger gcd = Gcd(n_, val.n_);
    val.n_ /= gcd;
    n_ /= gcd;
    gcd = Gcd(m_, val.m_);
    m_ /= gcd;
    val.m_ /= gcd;
    n_ *= val.m_;
    m_ *= val.n_;
    FixRational();
    return *this;
  }

  std::string toString() const {
    std::string result;
    result += n_.toString();
    if (m_ != 1) {
      result += '/';
      result += m_.toString();
    }
    return result;
  }


  std::string asDecimal(size_t precision = 0) const {
    if (precision == 0) {
      return (n_ / m_).toString();
    }
    BigInteger big_ans = n_;
    if (!big_ans.GetSign()) {
      big_ans = -big_ans;
    }
    std::string ans;
    std::string temp_val;
    big_ans *= BigInteger::BinPow(10, precision);
    big_ans /= m_;
    ans = big_ans.toString();
    if (ans.size() <= precision) {
      for (int i = 0; i < int(precision) - int(ans.size()); ++i) {
        temp_val += "0";
      }
      ans = "0." + temp_val + ans;
    } else {
      ans = ans.substr(0, ans.size() - precision) + "." + ans.substr(ans.size() - precision, ans.size() - 1);
    }
    if (n_ < 0) {
      ans = "-" + ans;
    }
    return ans;
  }

  explicit operator double() const {
    return static_cast<double>(n_) / static_cast<double>(m_);
  }

 private:
  BigInteger n_{0};
  BigInteger m_{1};

  void FixRational() {
    if (n_ == 0) {
      m_ = 1;
      return;
    }
    if (m_ < 0) {
      m_ = -m_;
      n_ = -n_;
    }
  }
};

bool operator==(const Rational &first, const Rational &second) {
  return first.GetM() == second.GetM() && first.GetN() == second.GetN();
}

bool operator!=(const Rational &first, const Rational &second) {
  return !(first == second);
}

bool operator<(const Rational &first, const Rational &second) {
  return first.GetN() * second.GetM() < first.GetM() * second.GetN();
}

bool operator>(const Rational &first, const Rational &second) {
  return first.GetN() * second.GetM() < first.GetM() * second.GetN();
}

bool operator<=(const Rational &first, const Rational &second) {
  return !(first > second);
}

bool operator>=(const Rational &first, const Rational &second) {
  return !(first < second);
}

Rational operator+(const Rational &first, const Rational &second) {
  Rational temp = first;
  temp += second;
  return temp;
}

Rational operator-(const Rational &first, const Rational &second) {
  Rational temp = first;
  temp -= second;
  return temp;
}

Rational operator*(const Rational &first, const Rational &second) {
  Rational temp = first;
  temp *= second;
  return temp;
}

Rational operator/(const Rational &first, const Rational &second) {
  Rational temp = first;
  temp /= second;
  return temp;
}

std::ostream &operator<<(std::ostream &out, const Rational &val) {
  std::string num = val.toString();
  out << num;
  return out;
}

std::istream &operator>>(std::istream &in, Rational &num) {
  std::string temp;
  in >> temp;
  size_t len = temp.length();
  size_t div = temp.find('/');
  if (div == std::string::npos) {
    num.n_ = temp;
    num.m_ = 1;
  } else {
    num.n_ = temp.substr(0, div);
    num.m_ = temp.substr(div + 1, len - div - 1);
  }
  return in;
}


template<int N>
class Finite {
 public:
  Finite() = default;

  Finite(int val) : number_((N + (val % N)) % N) {}

  long long Get_Num() const {
    return number_;
  }

  Finite<N> operator-() const {
    return Finite(-number_);
  }

  Finite<N> &operator+=(const Finite<N> &second_num) {
    number_ += second_num.Get_Num();
    number_ %= N;
    return *this;
  }

  Finite<N> &operator-=(const Finite<N> &second_num) {
    *this += -second_num;
    return *this;
  }

  Finite<N> &operator*=(const Finite<N> &second_num) {
    number_ *= second_num.Get_Num();
    number_ %= N;
    return *this;
  }

  Finite<N> &operator/=(const Finite<N> &second_num) {
    Finite<N> reverse_num = BinPow(second_num, N - 1);
    *this *= reverse_num;
    return *this;
  }

  Finite<N> &operator++() {
    *this += 1;
    return *this;
  }

  Finite<N> operator++(int) {
    Finite<N> temp = *this;
    ++*this;
    return temp;
  }

  Finite<N> &operator--() {
    *this -= 1;
    return *this;
  }

  Finite<N> operator--(int) {
    Finite<N> temp = *this;
    --*this;
    return temp;
  }

 protected:
  int number_{0};
};

template<int N>
std::ostream &operator<<(std::ostream &out, Finite<N> &num) {
  int temp = num.Get_Num();
  out << temp;
  return out;
}

template<int N>
static Finite<N> BinPow(Finite<N> num, size_t degree) {
  Finite<N> res = 1;
  while (degree) {
    if (degree % 2 == 1) {
      res = res * num;
    }
    num = num * num;
    degree /= 2;
  }
  return res;
}

template<int N>
Finite<N> operator+(const Finite<N> num_first, const Finite<N> num_second) {
  return Finite<N>((num_first.Get_Num() + num_second.Get_Num()) % N);
}

template<int N>
Finite<N> operator-(const Finite<N> num_first, const Finite<N> num_second) {
  return Finite<N>((num_first.Get_Num() - num_second.Get_Num()) % N);
}

template<int N>
Finite<N> operator*(const Finite<N> num_first, const Finite<N> num_second) {
  return Finite<N>((num_first.Get_Num() * num_second.Get_Num()) % N);
}

template<int N>
Finite<N> operator/(const Finite<N> num_first, const Finite<N> num_second) {
  Finite<N> reverse_num = BinPow(num_second, N - 2);
  return Finite<N>((num_first.Get_Num() * reverse_num.Get_Num()) % N);
}

template<int N>
bool operator==(const Finite<N> num_first, const Finite<N> num_second) {
  return num_first.Get_Num() == num_second.Get_Num();
}

template<int N>
bool operator!=(const Finite<N> num_first, const Finite<N> num_second) {
  return !(num_first == num_second);
}


template<unsigned N, unsigned M, typename Field = Rational>
class Matrix {
 public:
  Matrix() {
    std::vector<std::vector<Field>> temp(N, std::vector<Field>(M, 0));
    for (size_t i = 0; i < std::min(N, M); ++i) {
      temp[i][i] = 1;
    }
    elements_ = temp;
  }

  Matrix(const std::vector<std::vector<Field>> vec) : elements_(vec) {}

  Matrix(const std::vector<std::vector<int>> vec) {
    elements_ = std::vector<std::vector<Field>>(N, std::vector<Field>(M));
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        elements_[i][j] = static_cast<Field>(vec[i][j]);
      }
    }
  }

  std::vector<Field> &operator[](size_t i) {
    return this->elements_[i];
  }

  std::vector<Field> operator[](size_t i) const {
    return this->elements_[i];
  }

  Matrix<N, M, Field> &operator=(const Matrix<N, M, Field> &second_mat) {
    elements_ = second_mat.Get_Elements();
    return *this;
  }

  Matrix<N, M, Field> &operator+=(const Matrix<N, M, Field> &second_mat) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        elements_[i][j] += second_mat[i][j];
      }
    }
    return *this;
  }

  Matrix<N, M, Field> &operator-=(const Matrix<N, M, Field> &second_mat) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        elements_[i][j] -= second_mat[i][j];
      }
    }
    return *this;
  }

  Matrix<N, M, Field> &operator*=(const Field val) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        elements_[i][j] *= val;
      }
    }
    return *this;
  }

  Matrix<N, M, Field> &operator/=(const Field val) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        elements_[i][j] /= val;
      }
    }
    return *this;
  }

  std::vector<std::vector<Field>> Get_Elements() const {
    return elements_;
  }

  Matrix<M, N, Field> transposed() const {
    Matrix<M, N, Field> ans;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        ans[j][i] = elements_[i][j];
      }
    }
    return ans;
  }

  size_t rank() const {
    size_t rank = 0;
    std::vector<std::vector<Field>> temp_vec = elements_;
    for (size_t i = 0; i < M; ++i) {
      size_t flag = N;
      for (size_t j = rank; j < N; ++j) {
        if (temp_vec[j][i] != static_cast<Field>(0)) {
          flag = j;
          break;
        }
      }
      if (flag == N) {
        continue;
      }
      if (flag != rank) {
        std::swap(temp_vec[rank], temp_vec[flag]);
      }
      ++rank;
      for (size_t j = rank; j < N; ++j) {
        Field coefficient = -temp_vec[j][i] / temp_vec[i][i];
        for (size_t r = 0; r < M; ++r) {
          temp_vec[j][r] += temp_vec[i][r] * coefficient;
        }
      }
    }
    return rank;

  }

  std::vector<Field> getRow(unsigned i) const {
    return elements_[i];
  }

  std::vector<Field> getColumn(unsigned j) const {
    std::vector<Field> ans;
    for (size_t i = 0; i < N; ++i) {
      ans.push_back(elements_[i][j]);
    }
    return ans;
  }


  Matrix<N, N, Field> &operator*=(const Matrix<N, N, Field> second_mat) {
    *this = *this * second_mat;
    return *this;
  }

  Field det() const {
    static_assert(M == N, "CE");
    size_t cnt_swaps = 0;
    Field det = static_cast<Field>(1);
    std::vector<std::vector<Field>> temp_vec = elements_;
    for (size_t i = 0; i < N; ++i) {
      size_t flag = N;
      for (size_t j = i; j < N; ++j) {
        if (temp_vec[j][i] != static_cast<Field>(0)) {
          flag = j;
          break;
        }
      }
      if (flag == N) {
        return 0;
      }
      if (flag != i) {
        ++cnt_swaps;
        std::swap(temp_vec[i], temp_vec[flag]);
      }
      det *= temp_vec[i][i];
      for (size_t j = i + 1; j < N; ++j) {
        Field coefficient = -temp_vec[j][i] / temp_vec[i][i];
        for (size_t r = i + 1; r < N; ++r) {
          temp_vec[j][r] += temp_vec[i][r] * coefficient;
        }
      }
    }
    if (cnt_swaps % 2 == 1) {
      det = -det;
    }
    return det;
  }

  Field trace() const {
    Field ans = static_cast<Field>(0);
    for (size_t i = 0; i < std::min(N, M); ++i) {
      ans += elements_[i][i];
    }
    return ans;
  }

  Matrix<N, N, Field> inverted() const {
    std::vector<std::vector<Field>> temp_vec(N, std::vector<Field>(2 * N, static_cast<Field>(0)));
    std::vector<std::vector<Field>> ans(N, std::vector<Field>(N));
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        temp_vec[i][j] = elements_[i][j];
      }
    }
    for (size_t j = 0; j < N; ++j) {
      temp_vec[j][j + N] = static_cast<Field>(1);
    }
    for (size_t i = 0; i < N; ++i) {
      size_t flag = i;
      for (size_t j = i; j < N; ++j) {
        if (temp_vec[j][i] != static_cast<Field>(0)) {
          flag = j;
          break;
        }
      }
      if (flag != i) {
        std::swap(temp_vec[i], temp_vec[flag]);
      }
      for (size_t j = 0; j < N; ++j) {
        Field coefficient = -temp_vec[j][i] / temp_vec[i][i];
        if (i != j) {
          for (size_t r = i; r < 2 * N; ++r) {
            temp_vec[j][r] += temp_vec[i][r] * coefficient;
          }
        }
      }
    }
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < N; ++j) {
        ans[i][j] = temp_vec[i][j + N] / temp_vec[i][i];
      }
    }
    return Matrix<N, N, Field>(ans);
  }

  void invert() {
    *this = inverted();
  }

 protected:
  std::vector<std::vector<Field>> elements_;
};

template<unsigned N, unsigned M, typename Field = Rational>
bool operator==(const Matrix<N, M, Field> first_mat, const Matrix<N, M, Field> second_mat) {
  return first_mat.Get_Elements() == second_mat.Get_Elements();
}

template<unsigned N, unsigned M, typename Field = Rational>
bool operator!=(const Matrix<N, M, Field> first_mat, const Matrix<N, M, Field> second_mat) {
  return !(first_mat == second_mat);
}

template<unsigned N, unsigned M, typename Field = Rational>
Matrix<N, M, Field> operator/(Matrix<N, M, Field> mat, const Field num) {
  mat /= num;
  return mat;
}

template<unsigned N, unsigned M, typename Field = Rational>
Matrix<N, M, Field> operator*(Matrix<N, M, Field> mat, const Field num) {
  mat *= num;
  return mat;
}

template<unsigned N, unsigned M, typename Field = Rational>
Matrix<N, M, Field> operator*(const Field num, Matrix<N, M, Field> mat) {
  mat *= num;
  return mat;
}

template<unsigned N, unsigned M, typename Field = Rational>
Matrix<N, M, Field> operator+(const Matrix<N, M, Field> first_mat, const Matrix<N, M, Field> second_mat) {
  Matrix<N, M, Field> temp = first_mat;
  temp += second_mat;
  return temp;
}

template<unsigned N, unsigned M, typename Field = Rational>
Matrix<N, M, Field> operator-(const Matrix<N, M, Field> first_mat, const Matrix<N, M, Field> second_mat) {
  Matrix<N, M, Field> temp = first_mat;
  temp -= second_mat;
  return temp;
}

template<unsigned N, unsigned K, unsigned M, typename Field = Rational>
Matrix<N, M, Field> operator*(const Matrix<N, K, Field> &first_mat, const Matrix<K, M, Field> &second_mat) {
  Matrix<N, M, Field> ans;
  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < M; ++j) {
      ans[i][j] = static_cast<Field>(0);
      for (size_t r = 0; r < K; ++r) {
        ans[i][j] += first_mat[i][r] * second_mat[r][j];
      }
    }
  }
  return ans;
}


template<unsigned N, typename Field = Rational>
class SquareMatrix : public Matrix<N, N, Field> {
 public:

  SquareMatrix() = default;

  SquareMatrix(std::vector<std::vector<Field>> &vec) : Matrix<N, N, Field>(vec) {}

  SquareMatrix(std::vector<std::vector<int>> &vec) : Matrix<N, N, Field>(vec) {}
};


