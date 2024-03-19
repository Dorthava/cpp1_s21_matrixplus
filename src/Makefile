CC = g++
FLAGS = -std=c++17 -g -Wall -Werror -Wextra
SOURCES = s21_matrix_oop.cc tests.cc

ifeq ($(shell uname), Linux)
	TESTFLAGS=-lgtest -lgtest_main -lrt -lstdc++ -pthread -fprofile-arcs -ftest-coverage
else
	TESTFLAGS=-lgtest -lgtest_main -lpthread -fprofile-arcs -ftest-coverage -lstdc++
endif

all: test

s21_matrix_oop.a:
	$(CC) $(FLAGS) -c s21_matrix_oop.cc
	ar -src lib_s21_matrix_oop.a s21_matrix_oop.o
	rm -f *.o

test: clean
	$(CC) $(FLAGS) $(SOURCES) -o test $(TESTFLAGS) -lgcov --coverage
	./test

gcov_report: test
	lcov -t "gcov_report" --ignore-errors mismatch -o report.info --no-external -c -d .
	genhtml -o report report.info
	open ./report/index.html

test_leaks: test
	leaks --atExit -- ./test

valgrind: test
	CK_FORK=no valgrind --vgdb=no --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./test
	make clean

to_style:
	cp ../materials/linters/.clang-format ./
	clang-format -i $(SOURCES) s21_matrix_oop.h
	rm -rf .clang-format

style:
	cp ../materials/linters/.clang-format .clang-format
	clang-format -style=Google -n $(SOURCES) s21_matrix_oop.h

clean: 
	rm -rf *.o s21_matrix_oop.a report *.info *.gcda *.gcno *.gcov *.gch *.out *.txt test test.dSYM .clang-format