#ifndef __GEO_H__
#define __GEO_H__

#define NOMINMAX

#include <limits>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <thread>

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

constexpr double G_PI = 3.14159265358979323846;
constexpr double G_E = 2.71828;

constexpr int G_LOG_TRACE = 5;
constexpr int G_LOG_DEBUG = 4;
constexpr int G_LOG_WARN = 3;
constexpr int G_LOG_ERROR = 2;
constexpr int G_LOG_NONE = 0;

namespace geo {
	
	G_DLL_EXPORT int loglevel();

	G_DLL_EXPORT void loglevel(int x);

	template <class T>
	T maxvalue() {
		return std::numeric_limits<T>::max();
	}

	template <class T>
	T smallvalue() {
		return std::numeric_limits<T>::lowest();
	}
	
	template <class T>
	T minvalue() {
		return std::numeric_limits<T>::min();
	}

	template <class T>
	inline T min(T a, T b) {
		return a > b ? b : a;
	}

	template <class T>
	inline T max(T a, T b) {
		return a < b ? b : a;
	}

	template <class T>
	inline T sq(T a) {
		return a * a;
	}

	template <class T>
	inline T abs(T x) {
		return x < 0 ? -x : x;
	}

	template <class T>
	T deg(T x) {
		return x * 180.0 / G_PI;
	}

	template <class T>
	T rad(T x) {
		return x * G_PI / 180.0;
	}


	/**
	 * Monitors and controls the operation of a running program.
	 * Provides access to status-tracking functions, and cancelation flags.
	 */
	class Monitor {
	protected:
		bool m_cancel;
		float m_start;
		float m_end;
		float m_lastStatus;
		Monitor* m_monitor;

	public:

		/**
		 * \brief Initialize a Monitor whose progress range is 0-1 (0-100%).
		 */
		Monitor() : Monitor(nullptr, 0, 1) {}

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
			m_lastStatus(start),
			m_monitor(monitor) {
		}

		/**
		 * \brief Return true if the cancel flag is set.
		 */
		virtual bool canceled() const {
			if (m_monitor) {
				return m_monitor->canceled();
			}
			else {
				return m_cancel;
			}
		}

		virtual void cancel() {
			if (m_monitor) {
				m_monitor->cancel();
			}
			else {
				m_cancel = true;
			}
		}

		virtual void setCanceled(bool cancel) {
			if(m_monitor) {
				m_monitor->setCanceled(cancel);
			} else {
				m_cancel = cancel;
			}
		}

		virtual void status(float status, const std::string& message = "") {
			if(status < 0)
				status = m_lastStatus;
			m_lastStatus = status;
			std::cout << std::setprecision(1) << std::fixed << message << " " << (m_start + status * (m_end - m_start)) * 100.0f << "%\n";
		}

		virtual void error(const std::string& err) {
			std::cerr << err << "\n";
		}

		virtual void exception(const std::exception* ex) {
			error(ex == nullptr ? "Unknown exception" : ex->what());
		}

		virtual ~Monitor() {}
	};

	/**
	 * \brief Return a default Monitor object if none is available elsewhere.
	 *
	 * \return A pointer to a global Monitor object.
	 */
	Monitor* getDefaultMonitor();

} // geo


#define g_log(x, y) { if(geo::loglevel() <= y) std::cerr << std::setprecision(12) << x << std::endl; }

#define g_trace(x) g_log("TRACE  : " << x, G_LOG_TRACE)
#define g_debug(x) g_log("DEBUG  : " << x, G_LOG_DEBUG)
#define g_warn(x)  g_log("WARNING: " << x, G_LOG_WARN)
#define g_error(x) g_log("ERROR  : " << x, G_LOG_ERROR)

#define g_argerr(x) {std::stringstream _ss; _ss << x; throw std::invalid_argument(_ss.str());}
#define g_implerr(x) {std::stringstream _ss; _ss << x; throw std::runtime_error(_ss.str());}
#define g_runerr(x) {std::stringstream _ss; _ss << x; throw std::runtime_error(_ss.str());}


#endif


