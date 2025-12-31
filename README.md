# Book code for Competitive Programming

This repository holds a collection of **reference implementations for Competitive Programming** written in modern **C++ (≥ C++17)**.

The goal of this project is to provide a **clean, reusable, robust, and testable library** of commonly used data structures and algorithm templates that can be used in programming contests such as **Codeforces**, **AtCoder**, and **ICPC-style competitions**.

The code is designed to be:
- Easy to reuse during contests
- Well-tested with automated tests
- Structured for long-term maintenance and learning

---

## Build & Test

This project uses **CMake** and **Catch2** for building and testing.

### Requirements
- CMake ≥ 3.20
- A C++ compiler with C++20 support
- Git (for fetching dependencies)

### Configure & Build
```bash
cmake -S . -B build
cmake --build build
```

### Run all tests
```bash
ctest --test-dir build
```

### Run a subset of tests
```bash
ctest --test-dir build -R <pattern>
```

You can also run test executables directly from the build directory:
```bash
./build/tests/<binary_name>
```

All test cases are written using **Catch2**, and tests are automatically discovered and registered with CTest.

---

## License

All code in this repository is owned by the author unless otherwise noted.

---

## Contributing

Contributions are welcome.

You may contribute by:
- Adding new algorithm or data structure implementations
- Improving existing code
- Adding or extending test coverage
- Fixing bugs or improving documentation

General steps:
1. Fork the repository
2. Create a feature branch
3. Commit your changes
4. Open a Pull Request

Please ensure that all tests pass before submitting a PR.

---

## Acknowledgments

This project is inspired by the competitive programming community and existing open-source CP libraries, as well as shared algorithm references and contest experiences.

Thanks to all contributors and the broader competitive programming community for sharing knowledge and ideas.
