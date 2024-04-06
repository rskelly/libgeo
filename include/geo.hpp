#ifndef __GEO_H__
#define __GEO_H__

#define NOMINMAX

#include <limits>
#include <exception>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <thread>
#include <mutex>

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

namespace dijital {
	
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


	class Monitor;

	static std::mutex __monitor_mtx;
	static Monitor* __monitor;

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

		Monitor() :
			m_cancel(false),
			m_start(0),
			m_end(0),
			m_lastStatus(0) {}

	public:

		static Monitor& get() {
			// TODO: Race condition.
			if(!__monitor) {
				std::lock_guard<std::mutex> lk(__monitor_mtx);
				if(!__monitor)
					__monitor = new Monitor();
			}
			return *__monitor;
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
		void init(float start, float end) {
			m_cancel = false;
			m_start = start;
			m_end = end;
			m_lastStatus = start;
		}

		/**
		 * \brief Return true if the cancel flag is set.
		 */
		bool canceled() const {
			return m_cancel;
		}

		void cancel() {
			m_cancel = true;
		}

		void setCanceled(bool cancel) {
			m_cancel = cancel;
		}

		void status(float status, const std::string& message = "") {
			if(status < 0)
				status = m_lastStatus;
			m_lastStatus = status;
			std::cout << std::setprecision(1) << std::fixed << message << " " << (m_start + status * (m_end - m_start)) * 100.0f << "%\n";
		}

		void error(const std::string& err) {
			std::cerr << err << "\n";
		}

		void exception(const std::exception* ex) {
			error(ex == nullptr ? "Unknown exception" : ex->what());
		}

		~Monitor() {}
	};

} // geo


#define g_log(x, y) { if(dijital::loglevel() <= y) std::cerr << std::setprecision(12) << x << std::endl; }

#define g_trace(x) g_log("TRACE  : " << x, G_LOG_TRACE)
#define g_debug(x) g_log("DEBUG  : " << x, G_LOG_DEBUG)
#define g_warn(x)  g_log("WARNING: " << x, G_LOG_WARN)
#define g_error(x) g_log("ERROR  : " << x, G_LOG_ERROR)

#define g_argerr(x) {std::stringstream _ss; _ss << x; throw std::invalid_argument(_ss.str());}
#define g_implerr(x) {std::stringstream _ss; _ss << x; throw std::runtime_error(_ss.str());}
#define g_runerr(x) {std::stringstream _ss; _ss << x; throw std::runtime_error(_ss.str());}


#endif


