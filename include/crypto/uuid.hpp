/*
 * Adapted from https://gist.github.com/ne-sachirou/882192
 * std::rand() can be replaced with other algorithms as Xorshift for better perfomance
 * Random seed must be initialized by user
 *
 *  Created on: Nov 22, 2017
 *      Author: rob
 */

#ifndef LIBGEO_INCLUDE_CRYPTO_UUID_HPP_
#define LIBGEO_INCLUDE_CRYPTO_UUID_HPP_

#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

#include <string>

namespace geo {
namespace crypto {

/**
 * Provides UUID strings.
 */
class UUID {
public:

	/**
	 * Return a unique identifier string.
	 */
	static std::string uuid(){
		boost::uuids::uuid uuid = boost::uuids::random_generator()();
		return boost::uuids::to_string(uuid);
	}
};

}
}


#endif /* LIBGEO_INCLUDE_CRYPTO_UUID_HPP_ */
