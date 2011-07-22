/* -*- c++ -*- */
/*
 * Copyright 2010 Free Software Foundation, Inc.
 * 
 * This file is part of GNU Radio
 * 
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <gr_ofdm_alamouti_tx_cc.h>
#include <gr_io_signature.h>
#include <stdexcept>
#include <iostream>
#include <string.h>

gr_ofdm_alamouti_tx_cc_sptr
gr_make_ofdm_alamouti_tx_cc(int fft_length)
{
  return gr_ofdm_alamouti_tx_cc_sptr(new gr_ofdm_alamouti_tx_cc(fft_length));
}

gr_ofdm_alamouti_tx_cc::gr_ofdm_alamouti_tx_cc(int fft_length)
  : gr_sync_block("ofdm_alamouti_tx_cc",
		  gr_make_io_signature(1, 1, sizeof(gr_complex)*fft_length),
		  gr_make_io_signature(2, 2, sizeof(gr_complex)*fft_length)),
    d_fft_length(fft_length)
{
  // CAUTION: We force this to make sure the Alamouti code has two items to work with
  // to produce its two output items. Any blocks feeding this block should know that.
  set_output_multiple(2);
}


gr_ofdm_alamouti_tx_cc::~gr_ofdm_alamouti_tx_cc()
{
}

int
gr_ofdm_alamouti_tx_cc::work (int noutput_items,
			      gr_vector_const_void_star &input_items,
			      gr_vector_void_star &output_items)
{
  const gr_complex *in = (const gr_complex *) input_items[0];

  gr_complex *out_sym0 = (gr_complex *) output_items[0];
  gr_complex *out_sym1 = (gr_complex *) output_items[1];

  int iptr = 0, optr = 0;
  while(iptr < noutput_items-1) {
    // copy [s0, s1]; we start with i pointing to s0
    memcpy(&out_sym0[optr * d_fft_length],
	   &in[iptr * d_fft_length],
	   d_fft_length * sizeof(gr_complex));

    iptr++; // increment pointer to s1
    memcpy(&out_sym1[optr * d_fft_length],
	   &in[iptr * d_fft_length],
	   d_fft_length * sizeof(gr_complex));
    optr++; // we've produced 1 output symbol
    
    // copy [-s1*, s0*]
    // Start at s1
    for(int j = 0; j < d_fft_length; j++) {
      out_sym0[optr*d_fft_length + j] = -conj(in[iptr*d_fft_length + j]);
    }

    iptr--; // go back to s0 to output s0*
    for(int j = 0; j < d_fft_length; j++) {
      out_sym1[optr*d_fft_length + j] = conj(in[iptr*d_fft_length + j]);
    }

    optr++;  // we've produced the second output symbol
    iptr+=2; // we have used two full input symbols   
  }

  return optr;
}
