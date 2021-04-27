#include <iostream>
#include <string>
#include <cstring>
#include <vector>

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
    return static_cast<double>(num_[0]);
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
    int for_next = 0;
    int i;
    for (i = 0; i < static_cast<int>(val.num_.size()); ++i) {
      if (i >= static_cast<int>(num_.size())) {
        num_.push_back(0);
      }
      num_[i] += val.num_[i] + for_next;
      for_next = 0;
      if (num_[i] >= base_) {
        num_[i] -= base_;
        for_next = 1;
      }
    }
    while (for_next) {
      if (static_cast<int>(num_.size()) <= i) {
        num_.push_back(0);
      }
      ++num_[i];
      for_next = 0;
      if (num_[i] >= base_) {
        num_[i] -= base_;
        for_next = 1;
      }
      ++i;
    }
    Fix();
    return *this;
  }


  BigInteger &operator*=(const BigInteger &val) {
    BigInteger temp;
    temp.num_.resize(val.num_.size() + num_.size(), 0);
    for (size_t i = 0; i < val.num_.size(); ++i) {
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

  BigInteger &operator-=(const BigInteger &);

  BigInteger &operator%=(const BigInteger &);

  void ShiftRight() {
    if (num_.empty()) {
      num_.push_back(0);
      return;
    }
    num_.push_back(num_[num_.size() - 1]);
    for (size_t i = num_.size() - 2; i > 0; --i) {
      num_[i] = num_[i - 1];
    }
    num_[0] = 0;
  }

  BigInteger &operator/=(const BigInteger &);


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
  bool first_sign = first_val.GetSign();
  bool second_sign = second_val.GetSign();
  std::vector<int> first_num = first_val.GetNum();
  std::vector<int> second_num = second_val.GetNum();
  return (first_sign == second_sign && first_num == second_num);
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
  if (*this < 0 && val < 0) {
    BigInteger temp = -val;
    temp -= -*this;
    return *this = temp;
  }
  if (*this < 0) {
    BigInteger temp = val;
    temp += -*this;
    return *this = -temp;
  }
  if (val < 0) {
    return (*this += -val);
  }
  if (*this < val) {
    BigInteger temp = val;
    temp -= *this;
    return *this = -temp;
  }
  int for_next = 0;
  int i;
  for (i = 0; i < static_cast<int>(val.num_.size()); ++i) {
    if (static_cast<int>(num_.size()) <= i) {
      num_.push_back(0);
    }
    num_[i] -= val.num_[i] + for_next;
    for_next = 0;
    if (num_[i] < 0) {
      num_[i] += base_;
      for_next = 1;
    }
  }
  while (for_next) {
    if (static_cast<int>(num_.size()) <= i) {
      num_.push_back(0);
    }
    --num_[i];
    for_next = 0;
    if (num_[i] < 0) {
      num_[i] += base_;
      for_next = 1;
    }
    ++i;
  }
  Fix();
  return *this;
}

BigInteger &BigInteger::operator/=(const BigInteger &val) {
  BigInteger div = val;
  if (div < 0) {
    div = -div;
  }
  BigInteger result;
  BigInteger current_num;
  result.num_.resize(num_.size(), 0);
  for (int i = static_cast<int>(num_.size()) - 1; i >= 0; --i) {
    current_num.ShiftRight();
    current_num.num_[0] = num_[i];
    current_num.Fix();
    int temp_res = 0;
    int left = 0;
    int right = 1'000'000'000;
    while (left <= right) {
      int m = (left + right) / 2;
      BigInteger temp = m * div;
      if (temp <= current_num) {
        temp_res = m;
        left = m + 1;
      } else {
        right = m - 1;
      }
    }
    result.num_[i] = temp_res;
    current_num = current_num - div * temp_res;
  }
  result.sign_ = sign_ == val.sign_;
  result.Fix();
  return *this = result;
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

BigInteger Gcd(const BigInteger &first, const BigInteger &second) {
  BigInteger first_copy = first;
  BigInteger second_copy = second;
  if (first_copy < 0) {
    first_copy = -first_copy;
  }
  if (second_copy < 0) {
    second_copy = -second_copy;
  }
  while (second_copy) {
    first_copy %= second_copy;
    std::swap(first_copy, second_copy);
  }
  return first_copy;
}


class Rational {
 public:

  Rational() = default;

  Rational(const BigInteger &val) : n_(val), m_(1) {}

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
    FixRational();
    return *this;
  }

  Rational &operator-=(const Rational &val) {
    n_ = n_ * val.m_ - m_ * val.n_;
    m_ *= val.m_;
    FixRational();
    return *this;
  }

  Rational &operator*=(const Rational &val) {
    n_ *= val.n_;
    m_ *= val.m_;
    FixRational();
    return *this;
  }

  Rational &operator/=(const Rational &val) {
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
    if (big_ans < 0) {
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
    BigInteger gcd = Gcd(n_, m_);
    m_ /= gcd;
    n_ /= gcd;
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

