#ifndef __GEO_H__
#define __GEO_H__

#include <limits>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>

#ifdef __GNUC__
#define DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define DEPRECATED __declspec(deprecated)
#else
#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED
#endif

#ifdef _MSC_VER
#define G_DLL_EXPORT __declspec(dllexport)
#else
#define G_DLL_EXPORT
#endif

#define G_PI 3.14159265358979323846
#define G_E 2.71828

#define G_DBL_MAX_POS (std::numeric_limits<double>::max())
#define G_DBL_MAX_NEG (std::numeric_limits<double>::lowest())
#define G_DBL_MIN_POS (std::numeric_limits<double>::min())

#define G_FLT_MAX_POS (std::numeric_limits<float>::max())
#define G_FLT_MAX_NEG (std::numeric_limits<float>::lowest())

#define g_min(a, b) ((a) > (b) ? (b) : (a))
#define g_max(a, b) ((a) < (b) ? (b) : (a))
#define g_sq(a) ((a) * (a))
#define g_abs(x) ((x) < 0 ? -(x) : (x))
#define g_deg(x) ((x) * 180.0 / G_PI)
#define g_rad(x) ((x) * G_PI / 180.0)

G_DLL_EXPORT extern int g__loglevel;

#define G_LOG_TRACE 5
#define G_LOG_DEBUG 4
#define G_LOG_WARN 3
#define G_LOG_ERROR 2
#define G_LOG_NONE 0

#define g_loglevel(x) {g__loglevel = x;}

#define g_log(x, y) { if(g__loglevel <= y) std::cerr << std::setprecision(12) << x << std::endl; }

#define g_trace(x) g_log("TRACE:   " << x, G_LOG_TRACE)
#define g_debug(x) g_log("DEBUG:   " << x, G_LOG_DEBUG)
#define g_warn(x)  g_log("WARNING: " << x, G_LOG_WARN)
#define g_error(x) g_log("ERROR:   " << x, G_LOG_ERROR)

#define g_argerr(x) {std::stringstream _ss; _ss << x; throw std::invalid_argument(_ss.str());}
#define g_implerr(x) {std::stringstream _ss; _ss << x; throw std::runtime_error(_ss.str());}
#define g_runerr(x) {std::stringstream _ss; _ss << x; throw std::runtime_error(_ss.str());}


namespace geo {

/**
 * Monitors and controls the operation of a running program.
 * Provides access to status-tracking functions, and cancelation flags.
 */
class Monitor {
private:
	bool m_cancel;
	float m_start;
	float m_end;
	Monitor* m_monitor;

public:

	/**
	 * \brief Initialize a Monitor whose progress range is 0-1 (0-100%).
	 */
	Monitor() :
		m_cancel(false),
		m_start(0), m_end(1),
		m_monitor(nullptr) {
	}

	/**
	 * \brief Initialize a Monitor whose progress range is {start} to {end}.
	 *
	 * The status is transformed by status = start + status * (end - start).
	 *
	 * \param monitor The Monitor to which status updates are passed.
	 * \param start The starting status.
	 * \param end The ending status.
	 */
	Monitor(Monitor* monitor, float start, float end) :
		m_cancel(false),
		m_start(start), m_end(end),
		m_monitor(monitor) {
	}

	/**
	 * \brief Return true if the cancel flag is set.
	 */
	bool canceled() const {
		if(m_monitor) {
			return m_monitor->canceled();
		} else {
			return m_cancel;
		}
	}

	void cancel() {
		if(m_monitor) {
			m_monitor->cancel();
		} else {
			m_cancel = true;
		}
	}

	void status(float status, const std::string& message = "") {
		std::cout << std::setprecision(1) << std::fixed << message << " " << (m_start + status * (m_end - m_start)) * 100.0f << "%\n";
	}

	void error(const std::string& err) {
		std::cerr << err << "\n";
	}

	void exception(const std::exception* ex) {
		error(ex == nullptr ? "Unknown exception" : ex->what());
	}

};

/**
 * \brief Return a default Monitor object if none is available elsewhere.
 *
 * \return A pointer to a global Monitor object.
 */
Monitor* getDefaultMonitor();

} // geo


#endif


