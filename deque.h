#include<vector>
#include<iostream>

template<bool IsConst, typename Type>
class Types {
};

template<typename Type>
class Types<false, Type> {
 public:
  using pointer = Type *;
  using reference = Type &;
};

template<typename Type>
class Types<true, Type> {
 public:
  using pointer = const Type *;
  using reference = const Type &;
};

template<typename T>
class Deque {
 private:
  struct position;
 public:
  Deque() = default;

  Deque(int n) : Deque() {
    try {
      for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
        push_back(T());
      }
    } catch (...) {
      for (auto vec: pointers_) {
        if(vec != nullptr) {}
        delete[] vec;
      }
    }
  }

  Deque(int n, const T &item) : Deque() {
    for (size_t i = 0; i < static_cast<size_t>(n); ++i) {
      push_back(item);
    }
  }

  Deque(const Deque<T> &another) : Deque() {
    for (auto it: another) {
      push_back(it);
    }
  }

  ~Deque() {
    for (auto vec: pointers_) {
      if (vec != nullptr) {}
      delete[] vec;
    }
  }

  Deque<T> &operator=(Deque<T> &another) {
    pointers_ = std::vector<T *>(2, nullptr);
    begin_ = position(1, 0);
    end_ = position(1, 0);
    size_ = 0;
    for (auto it: another) {
      push_back(it);
    }
    return *this;
  }

  size_t size() const {
    return size_;
  }

  T &operator[](const int &n) {
    size_t vect_num = begin_.vector_pos + n / size_of_bucket_;
    size_t in_vect_num = begin_.in_vector_pos + n % size_of_bucket_;
    if (in_vect_num >= size_of_bucket_) {
      in_vect_num -= size_of_bucket_;
      ++vect_num;
    }
    return pointers_[vect_num][in_vect_num];
  }

  T operator[](const int &n) const {
    size_t vect_num = begin_.vector_pos + n / size_of_bucket_;
    size_t in_vect_num = begin_.in_vector_pos + n % size_of_bucket_;
    if (in_vect_num >= size_of_bucket_) {
      in_vect_num -= size_of_bucket_;
      ++vect_num;
    }
    return pointers_[vect_num][in_vect_num];
  }

  T &at(const int &n) {
    if (static_cast<size_t>(n) >= size_) {
      throw std::out_of_range("");
    }
    return *this[n];
  }

  T at(const int &n) const {
    if (n >= size_) {
      throw std::out_of_range("");
    }
    return *this[n];
  }

  void push_back(const T &item) {
    auto copy = pointers_;
    try {
      if (pointers_[end_.vector_pos] == nullptr) {
        pointers_[end_.vector_pos] = reinterpret_cast<T *>(new int8_t[sizeof(T) * size_of_bucket_]);
      }
      new(pointers_[end_.vector_pos] + end_.in_vector_pos) T(item);
      end_.increase();
      if (end_.vector_pos == pointers_.size()) {
        double_size();
      }
      ++size_;
    }
    catch (...) {
      if (end_.vector_pos == 1) {
        end_.decrease();
        if (pointers_[end_.vector_pos] != nullptr) {
          delete[] pointers_[end_.vector_pos];
        }
      }
      pointers_ = copy;
      throw;
    }
  }

  void pop_back() {
    end_.decrease();
    --size_;
  }

  void push_front(const T &item) {
    try {
      if (begin_.vector_pos == 0 && begin_.in_vector_pos == 0) {
        double_size();
      }
      begin_.decrease();
      if (pointers_[begin_.vector_pos] == nullptr) {
        pointers_[begin_.vector_pos] = reinterpret_cast<T *>(new int8_t[sizeof(T) * size_of_bucket_]);
      }
      new(pointers_[begin_.vector_pos] + begin_.in_vector_pos) T(item);
      ++size_;
    }
    catch (...) {
      begin_.increase();
      if (begin_.in_vector_pos == 0) {
        divide_size();
      }
      throw;
    }
  }

  void pop_front() {
    begin_.increase();
    --size_;
  }

  operator const Deque<T>() {

  }

  template<bool IsConst>
  struct iterator_template {
    friend class Deque<T>;

   private:
    using deq_type = std::conditional<IsConst, const Deque<T> *, Deque<T> *>;
    typename Types<IsConst, T>::pointer pointer;
    Deque<T>::position pos;
    const Deque<T> *deq;
   public:
    iterator_template(typename Types<IsConst, T>::pointer ptr, position ps, const Deque<T> *dq) :
            pointer(ptr), pos(ps),
            deq(dq) {}

    operator iterator_template<true>() {
      auto ans = iterator_template<true>(pos);
      return ans;
    }

    iterator_template<IsConst> &operator=(iterator_template<IsConst> other) {
      pointer = other.pointer;
      pos = other.pos;
      return *this;
    }

    typename Types<IsConst, T>::reference operator*() const {
      return *pointer;
    }

    typename Types<IsConst, T>::pointer operator->() const {
      return pointer;
    }

    bool operator==(iterator_template<IsConst> it1) {
      return it1.get_pos() == get_pos();
    }

    bool operator!=(iterator_template<IsConst> it1) {
      return it1.get_pos() != get_pos();
    }

    iterator_template<IsConst> operator+(const int n) {
      iterator_template<IsConst> temp = *this;
      temp += n;
      return temp;
    }

    iterator_template<IsConst> operator-(const int n) {
      iterator_template<IsConst> temp = *this;
      temp -= n;
      return temp;
    }


    bool operator>=(iterator_template<IsConst> it2) {
      return get_pos() >= it2.get_pos();
    }

    bool operator<=(iterator_template<IsConst> it2) {
      return get_pos() <= it2.get_pos();
    }

    bool operator>(iterator_template<IsConst> it2) {
      return get_pos() > it2.get_pos();
    }

    bool operator<(iterator_template<IsConst> it1) {
      return get_pos() < it1.get_pos();
    }

    int operator-(iterator_template<IsConst> it1) {
      return get_pos() - it1.get_pos();
    }

    size_t get_pos() {
      return pos.vector_pos * size_of_bucket_ + pos.in_vector_pos;
    }

    iterator_template<IsConst> &operator++() {
      pos.increase();
      if (pos.in_vector_pos != 0) {
        ++pointer;
      } else {
        pointer = deq->get_ptr(pos);
      }
      return *this;
    }

    iterator_template<IsConst> operator++(int) {
      auto it = *this;
      ++(*this);
      return it;
    }

    iterator_template<IsConst> &operator--() {
      pos.decrease();
      if (pos.in_vector_pos != size_of_bucket_ - 1) {
        --pointer;
      } else {
        pointer = deq->get_ptr(pos);
      }
      return *this;
    }

    iterator_template<IsConst> operator--(int) {
      return --*this;
    }

    iterator_template<IsConst> &operator+=(int n) {
      int temp = n + pos.vector_pos * size_of_bucket_ + pos.in_vector_pos;
      pos = position(temp);
      pointer = deq->get_ptr(pos);
      return *this;
    }

    iterator_template<IsConst> &operator-=(int n) {
      return (*this += -n);
    }
  };

  T *get_ptr(position pos) const {
    return &pointers_[pos.vector_pos][pos.in_vector_pos];
  }

  using const_iterator = iterator_template<true>;
  using iterator = iterator_template<false>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  iterator begin() {
    auto ptr = &pointers_[begin_.vector_pos][begin_.in_vector_pos];
    return iterator(ptr, begin_, this);
  }

  const_iterator begin() const {
    auto ptr = &pointers_[begin_.vector_pos][begin_.in_vector_pos];
    return const_iterator(ptr, begin_, this);
  }

  iterator end() {
    auto ptr = &pointers_[end_.vector_pos][end_.in_vector_pos];
    return iterator(ptr, end_, this);
  }

  const_iterator end() const {
    auto ptr = &pointers_[end_.vector_pos][end_.in_vector_pos];
    return const_iterator(ptr, end_, this);
  }

  const_iterator cbegin() const {
    return begin();
  }

  const_iterator cend() const {
    return end();
  }

  reverse_iterator rbegin() {
    position ps = end();
    ps.decrease();
    auto ptr = &pointers_[ps.vector_pos][ps.in_vector_pos];
    return reverse_iterator(ptr, ps, this);
  }

  const_reverse_iterator rbegin() const {
    position ps = end();
    ps.decrease();
    auto ptr = &pointers_[ps.vector_pos][ps.in_vector_pos];
    return const_reverse_iterator(ptr, ps, this);
  }

  reverse_iterator rend() {
    position ps = begin();
    ps.decrease();
    auto ptr = &pointers_[ps.vector_pos][ps.in_vector_pos];
    return reverse_iterator(ptr, ps);
  }

  const_reverse_iterator rend() const {
    position ps = begin();
    ps.decrease();
    auto ptr = &pointers_[ps.vector_pos][ps.in_vector_pos];
    return const_reverse_iterator(ptr, ps);
  }

  const_reverse_iterator crbegin() const {
    return rbegin();
  }

  const_reverse_iterator crend() const {
    return rend();
  }

  void insert(iterator given_iter, const T &item) {
    auto iter = end();
    try {
      --iter;
      push_back(*iter);
      for (; iter != given_iter; --iter) {
        *iter = *(iter - 1);
      }
      *iter = item;
    } catch (...) {
      ++iter;
      for (; iter != end(); ++iter) {
        *iter = *(iter + 1);
      }
      pop_back();
      throw;
    }
  }

  void erase(iterator given_iter) {
    auto reserve = *this;
    try {
      auto iter = given_iter;
      ++iter;
      for (; iter != end(); ++iter) {
        *(iter - 1) = *iter;
      }
      pop_back();
    } catch (...) {
      *this = reserve;
      throw;
    }
  }

 private:

  struct position {
    size_t vector_pos;
    size_t in_vector_pos;

    position() = default;

    position(size_t n) {
      vector_pos = n / size_of_bucket_;
      in_vector_pos = n % size_of_bucket_;
    }

    position(size_t first, size_t second) {
      vector_pos = first;
      in_vector_pos = second;
    }

    void increase() {
      if (in_vector_pos < size_of_bucket_ - 1) {
        ++in_vector_pos;
      } else {
        in_vector_pos = 0;
        ++vector_pos;
      }
    }

    void decrease() {
      if (in_vector_pos != 0) {
        --in_vector_pos;
      } else {
        in_vector_pos = size_of_bucket_ - 1;
        --vector_pos;
      }
    }
  };

  size_t size_ = 0;
  static const size_t size_of_bucket_ = 8;

  std::vector<T *> pointers_ = std::vector<T *>(2, nullptr);

  position begin_ = position(1, 0);
  position end_ = position(1, 0);

  void double_size() {
    size_t current_size = pointers_.size();
    size_t start = current_size / 2;
    std::vector<T *> temp_vector = std::vector<T *>(2 * current_size, nullptr);
    for (size_t i = 0; i < current_size; ++i) {
      temp_vector[i + start] = pointers_[i];
    }
    pointers_ = temp_vector;
    begin_.vector_pos += start;
    end_.vector_pos += start;
  }

  void divide_size() {
    size_t current_size = pointers_.size();
    size_t start = current_size / 2;
    std::vector<T *> temp_vector = std::vector<T *>(current_size / 2, nullptr);
    for (size_t i = 0; i < current_size; ++i) {
      temp_vector[i] = pointers_[i + start];
    }
    pointers_ = temp_vector;
    begin_.vector_pos -= start;
    end_.vector_pos -= start;
  }
};

