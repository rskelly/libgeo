/*
 * uuid.hpp
 *
 *  Created on: Nov 22, 2017
 *      Author: rob
 */

#ifndef LIBGEO_INCLUDE_CRYPTO_UUID_HPP_
#define LIBGEO_INCLUDE_CRYPTO_UUID_HPP_

#include <string>
#include <cstdlib>

//*Adapted from https://gist.github.com/ne-sachirou/882192
//*std::rand() can be replaced with other algorithms as Xorshift for better perfomance
//*Random seed must be initialized by user

namespace geo {
namespace crypto {

const std::string CHARS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

class UUID {
public:
	static std::string uuid(){
	  std::string uuid = std::string(36,' ');
	  int rnd = 0;
	  int r = 0;

	  uuid[8] = '-';
	  uuid[13] = '-';
	  uuid[18] = '-';
	  uuid[23] = '-';

	  uuid[14] = '4';

	  for(int i=0;i<36;i++){
		if (i != 8 && i != 13 && i != 18 && i != 14 && i != 23) {
		  if (rnd <= 0x02) {
			  rnd = 0x2000000 + (std::rand() * 0x1000000) | 0;
		  }
		  rnd >>= 4;
		  uuid[i] = CHARS[(i == 19) ? ((rnd & 0xf) & 0x3) | 0x8 : rnd & 0xf];
		}
	  }
	  return uuid;
	}
};

}
}


#endif /* LIBGEO_INCLUDE_CRYPTO_UUID_HPP_ */
