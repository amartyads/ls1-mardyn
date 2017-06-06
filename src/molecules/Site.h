#ifndef SITE_H_
#define SITE_H_

#include "utils/mardyn_assert.h"

#include "utils/Logger.h"
#include "utils/xmlfileUnits.h"
#include <array>
#include <cmath>
#include <cstdint>

using Log::global_log;

/** Site
 *
 * Sites are the center of a interactions. Depending on the interaction type
 * there are specialised derived classes of this basic site class.
 *
 * @author Martin Bernreuther et al. (2009)
 */
class Site {
public:
	double rx() const { return _r[0]; }  /**< get x-coordinate of position vector */
	double ry() const { return _r[1]; }  /**< get y-coordinate of position vector */
	double rz() const { return _r[2]; }  /**< get z-coordinate of position vector */
	std::array<double, 3> r() const { return _r; }  /**< get position vector */
	double m() const { return _m; }         /**< get mass */

	/**
	 * set the d-th component of the position
	 */
	void setR(int d, double r) {
		mardyn_assert(d < 3);
		_r[d] = r;
	}

	void setM(double m) {
		_m = m;
	}

	virtual ~Site() {}

protected:
	/// Constructor
    Site(double x = 0., double y = 0., double z = 0., double m = 0.)
        : _m(m) {
            _r[0] = x;
            _r[1] = y;
            _r[2] = z;
        }

	std::array<double, 3> _r; /**< position coordinates */
	double _m;    /**< mass */
};


/** Lennard-Jones 12-6 center
 * @author Martin Bernreuther et al.
 *
 * Lennard-Jones 12-6 interaction site. The potential between two LJ centers of the same type is 
 * given by 
 *
 * \f$[
 *  U_\text{LJ} = \epsilon \left[ \left(\frac{r}{\sigma}\right)^{6} - \left(\frac{r}{\sigma}\right)^{12} \right]
 * \f]
 *
 * where $r$ is the distance between the two LJ centers. See potforce.h for the detailed implementation.
 */
class LJcenter : public Site {
public:
	LJcenter():
		Site(0., 0., 0., 1.), _eps(0.), _sigma(0.), _rc(1.), _uLJshift6(1.) {}
	/** Constructor
     *
     * \param[in] x        relative x coordinate
     * \param[in] y        relative y coordinate
     * \param[in] z        relative z coordinate
     * \param[in] m        mass    
     * \param[in] epsilon  interaction strength
     * \param[in] sigma    interaction diameter
     * \param[in] rc       cutoff radius
     * \param[in] shift    0. for full LJ potential
     */
	LJcenter(double x, double y, double z, double m, double eps, double sigma, double rc, double shift)
			: Site(x, y, z, m), _eps(eps), _sigma(sigma), _rc(rc), _uLJshift6(shift) { }

	void readXML(XMLfileUnits& xmlconfig) {
		xmlconfig.getNodeValueReduced("coords/x", _r[0]);
		xmlconfig.getNodeValueReduced("coords/y", _r[1]);
		xmlconfig.getNodeValueReduced("coords/z", _r[2]);
		xmlconfig.getNodeValueReduced("mass", _m);
		xmlconfig.getNodeValueReduced("epsilon", _eps);
		xmlconfig.getNodeValueReduced("sigma", _sigma);
		if( xmlconfig.getNodeValueReduced("cutoff", _rc) == 0 ) {
			/* TODO: remove rc from LJ site? */
			_rc = -1.; /* set to invalid value  */
			global_log->warning() << "Cutoff radius for LJ site not specified" << std::endl;
		}
		_uLJshift6 = 0.0;
		xmlconfig.getNodeValueReduced("shifted", _uLJshift6);
	}
	
	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _m << "\t" << _eps << " " << _sigma << " " << _rc << " " << _uLJshift6;
	}
    /* TODO rename to epsilon */
	double eps() const { return _eps; }     /**< get interaction strength */
	double sigma() const { return _sigma; } /**< get interaction diameter */

	void setEps(double eps) {
		_eps = eps;
	}

	void setSigma(double sigma) {
		_sigma = sigma;
	}

	void setRC(double rc) {
		_rc = rc;
	}

	void setULJShift6(double uLJshift6) {
		_uLJshift6 = uLJshift6;
	}

    /* TODO: The following method is never used */
	bool TRUNCATED_SHIFTED() { return (_uLJshift6 != 0.0); } /**< get truncation option */

	double shift6() const { return _uLJshift6; }             /**< get energy shift of interaction potential */

private:
	double _eps;    /**< interaction strength */
	double _sigma;  /**< interaction diameter */

	// cutoff radius
	// it seems to me as if this is not the cutoff-radius which is used for the linked cells,
	// but a molecule specific one to determine the cutoff correction
	// TODO why may they be different!!!???
	double _rc;     /**< cutoff radius */
	double _uLJshift6; /**< truncation offset of the LJ potential */
};

/** Charge center
 */
class Charge : public Site {
public:
	Charge() :
			Site(0., 0., 0., 1.), _q(0.) {
	}
    /** Constructor
     *
     * \param[in] x        relative x coordinate
     * \param[in] y        relative y coordinate
     * \param[in] z        relative z coordinate
     * \param[in] m        mass    
     * \param[in] q        charge
     */
    Charge(double x, double y, double z, double m, double q)
			: Site(x, y, z, m), _q(q) { }

	void readXML(XMLfileUnits& xmlconfig) {
		xmlconfig.getNodeValueReduced("coords/x", _r[0]);
		xmlconfig.getNodeValueReduced("coords/y", _r[1]);
		xmlconfig.getNodeValueReduced("coords/z", _r[2]);
		xmlconfig.getNodeValueReduced("mass", _m);
		xmlconfig.getNodeValueReduced("charge", _q);
	}
	
	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _m << " " << _q;
	}
	double q() const { return _q; }  /**< get charge */

	void setQ(double q) {
		_q = q;
	}

private:
	double _q;  /**< charge */
};


/** Oriented site
 * @author Martin Bernreuther
 */
class OrientedSite : public Site {
public:
	double ex() const { return _e[0]; }
	double ey() const { return _e[1]; }
	double ez() const { return _e[2]; }
	std::array<double, 3> e() const { return _e; }  /**< Get pointer to the normalized orientation vector. */

	/** set the d-th component of the orientation vector */
	void setE(int d, double e) {
		mardyn_assert(d < 3);
		_e[d] = e;
	}

protected:
	/// Constructor
	OrientedSite(double x = 0., double y = 0., double z = 0., double m = 0., double ex = 0., double ey = 0., double ez = 0.)
			: Site(x, y, z, m) {
		_e[0] = ex;
		_e[1] = ey;
		_e[2] = ez;
	}

	std::array<double, 3> _e;  /**< Normalized orientation vector */

	void setOrientationVectorByPolarAngles(const double& theta_deg, const double& phi_deg)
	{
		const double fac = M_PI / 180.;  // translate: angle --> rad
		double theta_rad = theta_deg * fac;
		double phi_rad   = phi_deg   * fac;
		_e[0] = sin(theta_rad) * cos(phi_rad);
		_e[1] = sin(theta_rad) * sin(phi_rad);
		_e[2] = cos(theta_rad);
	}
};


/** Dipole
 * @author Martin Bernreuther
 *
 */
class Dipole : public OrientedSite {
public:
	Dipole():OrientedSite(0., 0., 0., 0., 1., 1., 1.), _absMy(0.) {}
    /** Constructor
     *
     * \param[in] x        relative x coordinate
     * \param[in] y        relative y coordinate
     * \param[in] z        relative z coordinate
     * \param[in] eMyx     x coordinate of the dipole moments normal
     * \param[in] eMyy     y coordinate of the dipole moments normal
     * \param[in] eMyz     z coordinate of the dipole moments normal
     * \param[in] absQ     dipole moments absolute value
     */
    Dipole(double x, double y, double z, double eMyx, double eMyy, double eMyz, double absMy)
			: OrientedSite(x, y, z, 0., eMyx, eMyy, eMyz), _absMy(absMy) {
	}

	void readXML(XMLfileUnits& xmlconfig) {
		xmlconfig.getNodeValueReduced("coords/x", _r[0]);
		xmlconfig.getNodeValueReduced("coords/y", _r[1]);
		xmlconfig.getNodeValueReduced("coords/z", _r[2]);
		xmlconfig.getNodeValueReduced("dipolemoment/abs", _absMy);
		bool bAngleInput = true;
		double theta, phi;
		bAngleInput = bAngleInput && xmlconfig.getNodeValueReduced("dipolemoment/theta", theta);
		bAngleInput = bAngleInput && xmlconfig.getNodeValueReduced("dipolemoment/phi",   phi);
		if(true == bAngleInput)
			this->setOrientationVectorByPolarAngles(theta, phi);
		else
		{
			xmlconfig.getNodeValueReduced("dipolemoment/x", _e[0]);
			xmlconfig.getNodeValueReduced("dipolemoment/y", _e[1]);
			xmlconfig.getNodeValueReduced("dipolemoment/z", _e[2]);
		}
		/* TODO normalization check */
	}

	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _e[0] << " " << _e[1] << " " << _e[2] << "\t" << _absMy;
	}

	double absMy() const { return _absMy; }  /**< Get the absolute value of the dipole moment. */

	/** set the value of the dipole moment */
	void setAbyMy(double my) {
		_absMy = my;
	}

private:
    /* TODO: move abs to oriented site. */
	double _absMy;  /**< absolute value of the dipole moment. */
};

/** Quadrupole
 * @author Martin Bernreuther
 */
class Quadrupole : public OrientedSite {
public:
	Quadrupole():OrientedSite(0., 0., 0., 0., 0., 0., 0.), _absQ(0.) {}
    /** Constructor
     *
     * \param[in] x        relative x coordinate
     * \param[in] y        relative y coordinate
     * \param[in] z        relative z coordinate
     * \param[in] eQx      x coordinate of the quadrupole moments normal
     * \param[in] eQy      y coordinate of the quadrupole moments normal
     * \param[in] eQz      z coordinate of the quadrupole moments normal
     * \param[in] absQ     quadrupole moments absolute value
     */
	Quadrupole(double x, double y, double z, double eQx, double eQy, double eQz, double absQ)
			: OrientedSite(x, y, z, 0., eQx, eQy, eQz), _absQ(absQ) {
	}

	void readXML(XMLfileUnits& xmlconfig) {
		xmlconfig.getNodeValueReduced("coords/x", _r[0]);
		xmlconfig.getNodeValueReduced("coords/y", _r[1]);
		xmlconfig.getNodeValueReduced("coords/z", _r[2]);
		xmlconfig.getNodeValueReduced("quadrupolemoment/abs", _absQ);
		bool bAngleInput = true;
		double theta, phi;
		bAngleInput = bAngleInput && xmlconfig.getNodeValueReduced("quadrupolemoment/theta", theta);
		bAngleInput = bAngleInput && xmlconfig.getNodeValueReduced("quadrupolemoment/phi",   phi);
		if(true == bAngleInput)
			this->setOrientationVectorByPolarAngles(theta, phi);
		else
		{
			xmlconfig.getNodeValueReduced("quadrupolemoment/x", _e[0]);
			xmlconfig.getNodeValueReduced("quadrupolemoment/y", _e[1]);
			xmlconfig.getNodeValueReduced("quadrupolemoment/z", _e[2]);
		}
		/* TODO normalization check */
	}
	
	/// write to stream
	void write(std::ostream& ostrm) const {
		ostrm << _r[0] << " " << _r[1] << " " << _r[2] << "\t" << _e[0] << " " << _e[1] << " " << _e[2] << " " << _absQ;
	}
	double absQ() const { return _absQ; }  /**< Get the absolute value of the quadrupole moment. */

	/** set the absolute value of teh quadrupole moment */
	void setAbsQ(double q) {
		_absQ = q;
	}

private:
    /* TODO: move abs to oriented site. */
	double _absQ;  /**< absolute value of the quadrupole moment. */
};

#endif  /* SITE_H_ */
