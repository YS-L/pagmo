/*****************************************************************************
 *   Copyright (C) 2004-2013 The PaGMO development team,                     *
 *   Advanced Concepts Team (ACT), European Space Agency (ESA)               *
 *   http://apps.sourceforge.net/mediawiki/pagmo                             *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Developers  *
 *   http://apps.sourceforge.net/mediawiki/pagmo/index.php?title=Credits     *
 *   act@esa.int                                                             *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.               *
 *****************************************************************************/

#include <boost/math/constants/constants.hpp>
#include <boost/lexical_cast.hpp>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>


#include "../exceptions.h"
#include "../types.h"
#include "base.h"
#include "cec2013.h"

#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625

namespace pagmo { namespace problem {

/// Constructor
/**
 * Will construct one of the 28 CEC2013 problems
 *
 * @param[in] fun_id The problem id. One of [1,2,...,28]
 * @param[in] d problem dimension. One of [2,5,10,20,30,...,100]
 * @param[in] dir The path where the CEC2013 input files are located.
 *                Two files are expected: "M_Dx.txt" and "shift_data.txt", where "x" is the problem dimension
 *
 * @see http://web.mysites.ntu.edu.sg/epnsugan/PublicSite/Shared%20Documents/CEC2013/cec13-c-code.zip to find
 * the files
 * @throws io_error if the files are not found
 */
cec2013::cec2013(unsigned int fun_id, problem::base::size_type d, const std::string& dir):base(d),m_problem_number(fun_id), m_y(d), m_z(d)
{
	if (!(d==2||d==5||d==10||d==20||d==30||d==40||d==50||d==60||d==70||d==80||d==90||d==100))
	{
		pagmo_throw(value_error, "Error: CEC2013 Test functions are only defined for dimensions 2,5,10,20,30,40,50,60,70,80,90,100.");
	}

	// We create the full file name for the rotation matrix
	std::string data_file_name(dir);
	data_file_name.append("M_D");
	data_file_name.append(boost::lexical_cast<std::string>(d));
	data_file_name.append(".txt");
	// And we read all datas into m_rotation_matrix
	{
	std::ifstream data_file(data_file_name.c_str());
	if (!data_file.is_open()) {
		pagmo_throw(io_error, std::string("Error: file not found. I was looking for (") + data_file_name.c_str() + ")");
	}
	std::istream_iterator<double> start(data_file), end;
	m_rotation_matrix = std::vector<double>(start,end);
	data_file.close();
	}

	// We create the full file name for the shift vector
	data_file_name = dir;
	data_file_name.append("shift_data.txt");
	// And we read all data into m_origin_shift
	{
	std::ifstream data_file(data_file_name.c_str());
	if (!data_file.is_open()) {
		pagmo_throw(io_error, std::string("Error: file not found. I was looking for ").append(data_file_name.c_str()));
	}
	std::istream_iterator<double> start(data_file), end;
	m_origin_shift = std::vector<double>(start,end);
	data_file.close();
	}
	// Set bounds. All CEC2013 problems have the same bounds
	set_bounds(-100,100);
}

/// Clone method.
base_ptr cec2013::clone() const
{
	return base_ptr(new cec2013(*this));
}

/// Implementation of the objective function.
void cec2013::objfun_impl(fitness_vector &f, const decision_vector &x) const
{
	size_type nx = get_dimension();
	switch(m_problem_number)
	{
	case 1:
		sphere_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],0);
		f[0]+=-1400.0;
		break;
	case 2:
		ellips_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=-1300.0;
		break;
	case 3:
		bent_cigar_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=-1200.0;
		break;
	case 4:
		discus_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=-1100.0;
		break;
	case 5:
		dif_powers_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],0);
		f[0]+=-1000.0;
		break;
	case 6:
		rosenbrock_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=-900.0;
		break;
	case 7:
		schaffer_F7_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=-800.0;
		break;
	case 8:
		ackley_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=-700.0;
		break;
	case 9:
		weierstrass_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=-600.0;
		break;
	case 10:
		griewank_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=-500.0;
		break;
	case 11:
		rastrigin_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],0);
		f[0]+=-400.0;
		break;
	case 12:
		rastrigin_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=-300.0;
		break;
	case 13:
		step_rastrigin_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=-200.0;
		break;
	case 14:
		schwefel_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],0);
		f[0]+=-100.0;
		break;
	case 15:
		schwefel_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=100.0;
		break;
	case 16:
		katsuura_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=200.0;
		break;
	case 17:
		bi_rastrigin_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],0);
		f[0]+=300.0;
		break;
	case 18:
		bi_rastrigin_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=400.0;
		break;
	case 19:
		grie_rosen_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=500.0;
		break;
	case 20:
		escaffer6_func(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=600.0;
		break;
	case 21:
		cf01(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=700.0;
		break;
	case 22:
		cf02(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],0);
		f[0]+=800.0;
		break;
	case 23:
		cf03(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=900.0;
		break;
	case 24:
		cf04(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=1000.0;
		break;
	case 25:
		cf05(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=1100.0;
		break;
	case 26:
		cf06(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=1200.0;
		break;
	case 27:
		cf07(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=1300.0;
		break;
	case 28:
		cf08(&x[0],&f[0],nx,&m_origin_shift[0],&m_rotation_matrix[0],1);
		f[0]+=1400.0;
		break;
	default:
		pagmo_throw(value_error, "Error: There are only 28 test functions in this test suite!");
		break;
	}

}

std::string cec2013::get_name() const
{
	std::string retval("CEC2013 - f");
	retval.append(boost::lexical_cast<std::string>(m_problem_number));
	switch(m_problem_number)
	{
	case 1:
		retval.append("(sphere_func)");
		break;
	case 2:
		retval.append("(ellips_func)");
		break;
	case 3:
		retval.append("(bent_cigar_func)");
		break;
	case 4:
		retval.append("(discus_func)");
		break;
	case 5:
		retval.append("(dif_powers_func_non_rotated)");
		break;
	case 6:
		retval.append("(rosenbrock_func)");
		break;
	case 7:
		retval.append("(schaffer_F7_func)");
		break;
	case 8:
		retval.append("(ackley_func)");
		break;
	case 9:
		retval.append("(weierstrass_func)");
		break;
	case 10:
		retval.append("(griewank_func)");
		break;
	case 11:
		retval.append("(rastrigin_func_non_rotated)");
		break;
	case 12:
		retval.append("(rastrigin_func)");
		break;
	case 13:
		retval.append("(step_rastrigin_func)");
		break;
	case 14:
		retval.append("(schwefel_func_non_rotated)");
		break;
	case 15:
		retval.append("(schwefel_func)");
		break;
	case 16:
		retval.append("(katsuura_func)");
		break;
	case 17:
		retval.append("(bi_rastrigin_func_non_rotated)");
		break;
	case 18:
		retval.append("(bi_rastrigin_func)");
		break;
	case 19:
		retval.append("(grie_rosen_func)");
		break;
	case 20:
		retval.append("(escaffer6_func)");
		break;
	case 21:
		retval.append("(cf01)");
		break;
	case 22:
		retval.append("(cf02)");
		break;
	case 23:
		retval.append("(cf03)");
		break;
	case 24:
		retval.append("(cf04)");
		break;
	case 25:
		retval.append("(cf05)");
		break;
	case 26:
		retval.append("(cf06)");
		break;
	case 27:
		retval.append("(cf07)");
		break;
	case 28:
		retval.append("(cf08)");
		break;
	default:
		pagmo_throw(value_error, "Error: There are only 28 test functions in this test suite!");
		break;
	}
	return retval;

}

void cec2013::sphere_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Sphere */
{
	shiftfunc(x, &m_y[0], nx, Os);
	if (r_flag==1)
		rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (int i=0; i<nx; i++)
		m_z[i]=m_y[i];
	f[0] = 0.0;
	for (int i=0; i<nx; i++)
	{
		f[0] += m_z[i]*m_z[i];
	}
}

void cec2013::ellips_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Ellipsoidal */
{
	int i;
	shiftfunc(x, &m_y[0], nx, Os);
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];
	oszfunc (&m_z[0], &m_y[0], nx);
	f[0] = 0.0;
	for (i=0; i<nx; i++)
	{
		f[0] += pow(10.0,6.0*i/(nx-1))*m_y[i]*m_y[i];
	}
}

void cec2013::bent_cigar_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Bent_Cigar */
{
	int i;
	double beta=0.5;
	shiftfunc(x, &m_y[0], nx, Os);
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];
	asyfunc (&m_z[0], &m_y[0], nx,beta);
	if (r_flag==1)
		rotatefunc(&m_y[0], &m_z[0], nx, &Mr[nx*nx]);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	f[0] = m_z[0]*m_z[0];
	for (i=1; i<nx; i++)
	{
		f[0] += pow(10.0,6.0)*m_z[i]*m_z[i];
	}
}

void cec2013::discus_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Discus */
{
	int i;
	shiftfunc(x, &m_y[0], nx, Os);
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];
	oszfunc (&m_z[0], &m_y[0], nx);

	f[0] = pow(10.0,6.0)*m_y[0]*m_y[0];
	for (i=1; i<nx; i++)
	{
		f[0] += m_y[i]*m_y[i];
	}
}

void cec2013::dif_powers_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Different Powers */
{
	int i;
	shiftfunc(x, &m_y[0], nx, Os);
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];
	f[0] = 0.0;
	for (i=0; i<nx; i++)
	{
		f[0] += pow(fabs(m_z[i]),2+4*i/(nx-1));
	}
	f[0]=pow(f[0],0.5);
}

void cec2013::rosenbrock_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Rosenbrock's */
{
	int i;
	double tmp1,tmp2;
	shiftfunc(x, &m_y[0], nx, Os);//shift
	for (i=0; i<nx; i++)//shrink to the orginal search range
	{
		m_y[i]=m_y[i]*2.048/100;
	}
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);//rotate
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];
	for (i=0; i<nx; i++)//shift to orgin
	{
		m_z[i]=m_z[i]+1;
	}

	f[0] = 0.0;
	for (i=0; i<nx-1; i++)
	{
		tmp1=m_z[i]*m_z[i]-m_z[i+1];
		tmp2=m_z[i]-1.0;
		f[0] += 100.0*tmp1*tmp1 +tmp2*tmp2;
	}
}

void cec2013::schaffer_F7_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Schwefel's 1.2  */
{
	int i;
	double tmp;
	shiftfunc(x, &m_y[0], nx, Os);
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];
	asyfunc (&m_z[0], &m_y[0], nx, 0.5);
	for (i=0; i<nx; i++)
		m_z[i] = m_y[i]*pow(10.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
	rotatefunc(&m_z[0], &m_y[0], nx, &Mr[nx*nx]);
	else
	for (i=0; i<nx; i++)
		m_y[i]=m_z[i];

	for (i=0; i<nx-1; i++)
		m_z[i]=pow(m_y[i]*m_y[i]+m_y[i+1]*m_y[i+1],0.5);
	f[0] = 0.0;
	for (i=0; i<nx-1; i++)
	{
	  tmp=sin(50.0*pow(m_z[i],0.2));
	  f[0] += pow(m_z[i],0.5)+pow(m_z[i],0.5)*tmp*tmp ;
	}
	f[0] = f[0]*f[0]/(nx-1)/(nx-1);
}

void cec2013::ackley_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Ackley's  */
{
	int i;
	double sum1, sum2;

	shiftfunc(x, &m_y[0], nx, Os);
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	asyfunc (&m_z[0], &m_y[0], nx, 0.5);
	for (i=0; i<nx; i++)
		m_z[i] = m_y[i]*pow(10.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
	rotatefunc(&m_z[0], &m_y[0], nx, &Mr[nx*nx]);
	else
	for (i=0; i<nx; i++)
		m_y[i]=m_z[i];

	sum1 = 0.0;
	sum2 = 0.0;
	for (i=0; i<nx; i++)
	{
		sum1 += m_y[i]*m_y[i];
		sum2 += cos(2.0*boost::math::constants::pi<double>()*m_y[i]);
	}
	sum1 = -0.2*sqrt(sum1/nx);
	sum2 /= nx;
	f[0] =  E - 20.0*exp(sum1) - exp(sum2) +20.0;
}

void cec2013::weierstrass_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Weierstrass's  */
{
	int i,j,k_max;
	double sum=0,sum2=0, a, b;

	shiftfunc(x, &m_y[0], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
	{
		m_y[i]=m_y[i]*0.5/100;
	}
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	asyfunc (&m_z[0], &m_y[0], nx, 0.5);
	for (i=0; i<nx; i++)
		m_z[i] = m_y[i]*pow(10.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
	rotatefunc(&m_z[0], &m_y[0], nx, &Mr[nx*nx]);
	else
	for (i=0; i<nx; i++)
		m_y[i]=m_z[i];

	a = 0.5;
	b = 3.0;
	k_max = 20;
	f[0] = 0.0;
	for (i=0; i<nx; i++)
	{
		sum = 0.0;
		sum2 = 0.0;
		for (j=0; j<=k_max; j++)
		{
			sum += pow(a,j)*cos(2.0*boost::math::constants::pi<double>()*pow(b,j)*(m_y[i]+0.5));
			sum2 += pow(a,j)*cos(2.0*boost::math::constants::pi<double>()*pow(b,j)*0.5);
		}
		f[0] += sum;
	}
	f[0] -= nx*sum2;
}

void cec2013::griewank_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Griewank's  */
{
	int i;
	double s, p;

	shiftfunc(x, &m_y[0], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
	{
		m_y[i]=m_y[i]*600.0/100.0;
	}
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	for (i=0; i<nx; i++)
		m_z[i] = m_z[i]*pow(100.0,1.0*i/(nx-1)/2.0);


	s = 0.0;
	p = 1.0;
	for (i=0; i<nx; i++)
	{
		s += m_z[i]*m_z[i];
		p *= cos(m_z[i]/sqrt(1.0+i));
	}
	f[0] = 1.0 + s/4000.0 - p;
}

void cec2013::rastrigin_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Rastrigin's  */
{
	int i;
	double alpha=10.0,beta=0.2;
	shiftfunc(x, &m_y[0], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
	{
		m_y[i]=m_y[i]*5.12/100;
	}

	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	oszfunc (&m_z[0], &m_y[0], nx);
	asyfunc (&m_y[0], &m_z[0], nx, beta);

	if (r_flag==1)
	rotatefunc(&m_z[0], &m_y[0], nx, &Mr[nx*nx]);
	else
	for (i=0; i<nx; i++)
		m_y[i]=m_z[i];

	for (i=0; i<nx; i++)
	{
		m_y[i]*=pow(alpha,1.0*i/(nx-1)/2);
	}

	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	f[0] = 0.0;
	for (i=0; i<nx; i++)
	{
		f[0] += (m_z[i]*m_z[i] - 10.0*cos(2.0*boost::math::constants::pi<double>()*m_z[i]) + 10.0);
	}
}

void cec2013::step_rastrigin_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Noncontinuous Rastrigin's  */
{
	int i;
	double alpha=10.0,beta=0.2;
	shiftfunc(x, &m_y[0], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
	{
		m_y[i]=m_y[i]*5.12/100;
	}

	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	for (i=0; i<nx; i++)
	{
		if (fabs(m_z[i])>0.5)
		m_z[i]=floor(2*m_z[i]+0.5)/2;
	}

	oszfunc (&m_z[0], &m_y[0], nx);
	asyfunc (&m_y[0], &m_z[0], nx, beta);

	if (r_flag==1)
	rotatefunc(&m_z[0], &m_y[0], nx, &Mr[nx*nx]);
	else
	for (i=0; i<nx; i++)
		m_y[i]=m_z[i];

	for (i=0; i<nx; i++)
	{
		m_y[i]*=pow(alpha,1.0*i/(nx-1)/2);
	}

	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	f[0] = 0.0;
	for (i=0; i<nx; i++)
	{
		f[0] += (m_z[i]*m_z[i] - 10.0*cos(2.0*boost::math::constants::pi<double>()*m_z[i]) + 10.0);
	}
}

void cec2013::schwefel_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Schwefel's  */
{
	int i;
	double tmp;
	shiftfunc(x, &m_y[0], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
	{
		m_y[i]*=1000/100;
	}
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	for (i=0; i<nx; i++)
		m_y[i] = m_z[i]*pow(10.0,1.0*i/(nx-1)/2.0);

	for (i=0; i<nx; i++)
		m_z[i] = m_y[i]+4.209687462275036e+002;

	f[0]=0;
	for (i=0; i<nx; i++)
	{
		if (m_z[i]>500)
		{
			f[0]-=(500.0-fmod(m_z[i],500))*sin(pow(500.0-fmod(m_z[i],500),0.5));
			tmp=(m_z[i]-500.0)/100;
			f[0]+= tmp*tmp/nx;
		}
		else if (m_z[i]<-500)
		{
			f[0]-=(-500.0+fmod(fabs(m_z[i]),500))*sin(pow(500.0-fmod(fabs(m_z[i]),500),0.5));
			tmp=(m_z[i]+500.0)/100;
			f[0]+= tmp*tmp/nx;
		}
		else
			f[0]-=m_z[i]*sin(pow(fabs(m_z[i]),0.5));
	}
	f[0]=4.189828872724338e+002*nx+f[0];
}

void cec2013::katsuura_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Katsuura  */
{
	int i,j;
	double temp,tmp1,tmp2,tmp3;
	tmp3=pow(1.0*nx,1.2);
	shiftfunc(x, &m_y[0], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
	{
		m_y[i]*=5.0/100.0;
	}
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	for (i=0; i<nx; i++)
		m_z[i] *=pow(100.0,1.0*i/(nx-1)/2.0);

	if (r_flag==1)
	rotatefunc(&m_z[0], &m_y[0], nx, &Mr[nx*nx]);
	else
	for (i=0; i<nx; i++)
		m_y[i]=m_z[i];

	f[0]=1.0;
	for (i=0; i<nx; i++)
	{
		temp=0.0;
		for (j=1; j<=32; j++)
		{
			tmp1=pow(2.0,j);
			tmp2=tmp1*m_y[i];
			temp += fabs(tmp2-floor(tmp2+0.5))/tmp1;
		}
		f[0] *= pow(1.0+(i+1)*temp,10.0/tmp3);
	}
	tmp1=10.0/nx/nx;
	f[0]=f[0]*tmp1-tmp1;

}

void cec2013::bi_rastrigin_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Lunacek Bi_rastrigin Function */
{
	int i;
	double mu0=2.5,d=1.0,s,mu1,tmp,tmp1,tmp2;
	double *tmpx;
	tmpx=(double *)malloc(sizeof(double)  *  nx);
	s=1.0-1.0/(2.0*pow(nx+20.0,0.5)-8.2);
	mu1=-pow((mu0*mu0-d)/s,0.5);

	shiftfunc(x, &m_y[0], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
	{
		m_y[i]*=10.0/100.0;
	}

	for (i = 0; i < nx; i++)
	{
		tmpx[i]=2*m_y[i];
		if (Os[i] < 0.)
			tmpx[i] *= -1.;
	}

	for (i=0; i<nx; i++)
	{
		m_z[i]=tmpx[i];
		tmpx[i] += mu0;
	}
	if (r_flag==1)
		rotatefunc(&m_z[0], &m_y[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_y[i]=m_z[i];

	for (i=0; i<nx; i++)
		m_y[i] *=pow(100.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
		rotatefunc(&m_y[0], &m_z[0], nx, &Mr[nx*nx]);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	tmp1=0.0;tmp2=0.0;
	for (i=0; i<nx; i++)
	{
		tmp = tmpx[i]-mu0;
		tmp1 += tmp*tmp;
		tmp = tmpx[i]-mu1;
		tmp2 += tmp*tmp;
	}
	tmp2 *= s;
	tmp2 += d*nx;
	tmp=0;
	for (i=0; i<nx; i++)
	{
		tmp+=cos(2.0*boost::math::constants::pi<double>()*m_z[i]);
	}

	if(tmp1<tmp2)
		f[0] = tmp1;
	else
		f[0] = tmp2;
	f[0] += 10.0*(nx-tmp);
	free(tmpx);
}

void cec2013::grie_rosen_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Griewank-Rosenbrock  */
{
	int i;
	double temp,tmp1,tmp2;

	shiftfunc(x, &m_y[0], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
	{
		m_y[i]=m_y[i]*5/100;
	}
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	for (i=0; i<nx; i++)//shift to orgin
	{
		m_z[i]=m_y[i]+1;
	}

	f[0]=0.0;
	for (i=0; i<nx-1; i++)
	{
		tmp1 = m_z[i]*m_z[i]-m_z[i+1];
		tmp2 = m_z[i]-1.0;
		temp = 100.0*tmp1*tmp1 + tmp2*tmp2;
		 f[0] += (temp*temp)/4000.0 - cos(temp) + 1.0;
	}
	tmp1 = m_z[nx-1]*m_z[nx-1]-m_z[0];
	tmp2 = m_z[nx-1]-1.0;
	temp = 100.0*tmp1*tmp1 + tmp2*tmp2;;
	 f[0] += (temp*temp)/4000.0 - cos(temp) + 1.0 ;
}

void cec2013::escaffer6_func (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Expanded Scaffer¡¯s F6  */
{
	int i;
	double temp1, temp2;
	shiftfunc(x, &m_y[0], nx, Os);
	if (r_flag==1)
	rotatefunc(&m_y[0], &m_z[0], nx, Mr);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	asyfunc (&m_z[0], &m_y[0], nx, 0.5);
	if (r_flag==1)
		rotatefunc(&m_y[0], &m_z[0], nx, &Mr[nx*nx]);
	else
	for (i=0; i<nx; i++)
		m_z[i]=m_y[i];

	f[0] = 0.0;
	for (i=0; i<nx-1; i++)
	{
		temp1 = sin(sqrt(m_z[i]*m_z[i]+m_z[i+1]*m_z[i+1]));
		temp1 =temp1*temp1;
		temp2 = 1.0 + 0.001*(m_z[i]*m_z[i]+m_z[i+1]*m_z[i+1]);
		f[0] += 0.5 + (temp1-0.5)/(temp2*temp2);
	}
	temp1 = sin(sqrt(m_z[nx-1]*m_z[nx-1]+m_z[0]*m_z[0]));
	temp1 =temp1*temp1;
	temp2 = 1.0 + 0.001*(m_z[nx-1]*m_z[nx-1]+m_z[0]*m_z[0]);
	f[0] += 0.5 + (temp1-0.5)/(temp2*temp2);
}

void cec2013::cf01 (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Composition Function 1 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10, 20, 30, 40, 50};
	double bias[5] = {0, 100, 200, 300, 400};

	i=0;
	rosenbrock_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/1e+4;
	i=1;
	dif_powers_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/1e+10;
	i=2;
	bent_cigar_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/1e+30;
	i=3;
	discus_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/1e+10;
	i=4;
	sphere_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],0);
	fit[i]=10000*fit[i]/1e+5;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cec2013::cf02 (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Composition Function 2 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {20,20,20};
	double bias[3] = {0, 100, 200};
	for(i=0;i<cf_num;i++)
	{
		schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	}
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cec2013::cf03 (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Composition Function 3 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {20,20,20};
	double bias[3] = {0, 100, 200};
	for(i=0;i<cf_num;i++)
	{
		schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	}
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cec2013::cf04 (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Composition Function 4 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {20,20,20};
	double bias[3] = {0, 100, 200};
	i=0;
	schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/4e+3;
	i=1;
	rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/1e+3;
	i=2;
	weierstrass_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/400;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cec2013::cf05 (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Composition Function 4 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {10,30,50};
	double bias[3] = {0, 100, 200};
	i=0;
	schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/4e+3;
	i=1;
	rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/1e+3;
	i=2;
	weierstrass_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/400;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cec2013::cf06 (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Composition Function 6 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,10,10,10,10};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/4e+3;
	i=1;
	rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/1e+3;
	i=2;
	ellips_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/1e+10;
	i=3;
	weierstrass_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/400;
	i=4;
	griewank_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=1000*fit[i]/100;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);

}

void cec2013::cf07 (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Composition Function 7 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,10,10,20,20};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	griewank_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/100;
	i=1;
	rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/1e+3;
	i=2;
	schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/4e+3;
	i=3;
	weierstrass_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/400;
	i=4;
	sphere_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],0);
	fit[i]=10000*fit[i]/1e+5;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cec2013::cf08 (const double *x, double *f, int nx, const double *Os,const double *Mr,int r_flag) const /* Composition Function 8 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,20,30,40,50};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	grie_rosen_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/4e+3;
	i=1;
	schaffer_F7_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/4e+6;
	i=2;
	schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/4e+3;
	i=3;
	escaffer6_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag);
	fit[i]=10000*fit[i]/2e+7;
	i=4;
	sphere_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],0);
	fit[i]=10000*fit[i]/1e+5;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cec2013::shiftfunc (const double *x, double *xshift, int nx, const double *Os) const {
	int i;
	for (i=0; i<nx; i++)
	{
		xshift[i]=x[i]-Os[i];
	}
}

void cec2013::rotatefunc (const double *x, double *xrot, int nx,const double *Mr) const
{
	int i,j;
	for (i=0; i<nx; i++)
	{
		xrot[i]=0;
			for (j=0; j<nx; j++)
			{
				xrot[i]=xrot[i]+x[j]*Mr[i*nx+j];
			}
	}
}

void cec2013::asyfunc (const double *x, double *xasy, int nx, double beta) const
{
	int i;
	for (i=0; i<nx; i++)
	{
		if (x[i]>0)
		xasy[i]=pow(x[i],1.0+beta*i/(nx-1)*pow(x[i],0.5));
	}
}

void cec2013::oszfunc (const double *x, double *xosz, int nx) const
{
	int i,sx;
	double c1,c2,xx=0;
	for (i=0; i<nx; i++)
	{
		if (i==0||i==nx-1)
		{
			if (x[i]!=0)
				xx=log(fabs(x[i]));
			if (x[i]>0)
			{
				c1=10;
				c2=7.9;
			}
			else
			{
				c1=5.5;
				c2=3.1;
			}
			if (x[i]>0)
				sx=1;
			else if (x[i]==0)
				sx=0;
			else
				sx=-1;
			xosz[i]=sx*exp(xx+0.049*(sin(c1*xx)+sin(c2*xx)));
		}
		else
			xosz[i]=x[i];
	}
}

void cec2013::cf_cal(const double *x, double *f, int nx, const double *Os,double * delta,double * bias,double * fit, int cf_num) const
{
	int i,j;
	double *w;
	double w_max=0,w_sum=0;
	w=(double *)malloc(cf_num * sizeof(double));
	for (i=0; i<cf_num; i++)
	{
		fit[i]+=bias[i];
		w[i]=0;
		for (j=0; j<nx; j++)
		{
			w[i]+=pow(x[j]-Os[i*nx+j],2.0);
		}
		if (w[i]!=0)
			w[i]=pow(1.0/w[i],0.5)*exp(-w[i]/2.0/nx/pow(delta[i],2.0));
		else
			w[i]=INF;
		if (w[i]>w_max)
			w_max=w[i];
	}

	for (i=0; i<cf_num; i++)
	{
		w_sum=w_sum+w[i];
	}
	if(w_max==0)
	{
		for (i=0; i<cf_num; i++)
			w[i]=1;
		w_sum=cf_num;
	}
	f[0] = 0.0;
	for (i=0; i<cf_num; i++)
	{
		f[0]=f[0]+w[i]/w_sum*fit[i];
	}
	free(w);
}

}} //namespaces

#undef EPS
#undef E
#undef INF

BOOST_CLASS_EXPORT_IMPLEMENT(pagmo::problem::cec2013);
