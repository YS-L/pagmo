#ifndef THROTTLE_H
#define THROTTLE_H

#include <numeric>

#include "../../serialization.h"
#include "../astro_constants.h"
#include "../epoch.h"

namespace kep_toolbox { namespace sims_flanagan{

/// A single throttle
/**
 * This class models a single throttle in the Sims-Flanagan model. It essentialy contains the cartesian
 * components of one throttle (non dimensional impulse)
 *impulse
 * @author David di Lorenzo
 * @author Dario Izzo (dario.izzo _AT_ googlemail.com)
 */

class throttle {
public:
	throttle():m_start(),m_end() {
		m_value[0] = 0;
		m_value[1] = 0;
		m_value[2] = 0;
	}

	throttle(epoch _start, epoch _end, const array3D& _value)
		: m_start(_start), m_end(_end), m_value(_value) {}

	epoch get_start() const {
		return m_start;
	}

	epoch get_end() const {
		return m_end;
	}

	const array3D& get_value() const {
		return m_value;
	}

	double get_norm() const {
		return std::sqrt(std::inner_product(m_value.begin(), m_value.end(), m_value.begin(), 0.));
	}


private:
	friend class boost::serialization::access;
	template <class Archive>
	void serialize(Archive &ar, const unsigned int)
	{
		ar & m_start;
		ar & m_end;
		ar & m_value;
	}
	epoch m_start;
	epoch m_end;
	array3D m_value;
};

}} //namespaces

#endif
