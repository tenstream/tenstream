/*--------------------------------------------------------------------
 * $Id: twomaxrnd.c 2623 2011-12-23 10:52:38Z nina.crnivec $
 *
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

/* full solution of the multi-layer twostream equation solar + thermal */
/* with maximum-random overlap assumption for partial cloudiness */
/* written by Nina Crnivec */

#ifndef __twomaxrnd_h
#define __twomaxrnd_h

#if defined (__cplusplus)
extern "C" {
#endif

    int twostream_maxrand (double *dtau_c, double *omega0_c, double *g_c, //"_c" for cloudy region
            double *dtau_f, double *omega0_f, double *g_f, //"_f" for cloud-free region
            double *cf, int nlev, double S0, double mu0, double Ag,
            double Bg, double *B, int delta, int flagSolar, int flagThermal,
            double *Edir, double *Edn, double *Eup);

#if defined (__cplusplus)
}
#endif

#endif /* _twomaxrnd_h */
