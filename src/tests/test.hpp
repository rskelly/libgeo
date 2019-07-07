/*
 * test.hpp
 *
 *  Created on: Jul 6, 2019
 *      Author: rob
 */

#ifndef SRC_TESTS_TEST_HPP_
#define SRC_TESTS_TEST_HPP_


namespace geo {
namespace test {

class test {
public:
	virtual bool run(int argc, char** argv) = 0;
	virtual ~test() {}
};
}
}


#endif /* SRC_TESTS_TEST_HPP_ */
