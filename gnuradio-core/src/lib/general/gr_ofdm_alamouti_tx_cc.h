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

#ifndef INCLUDED_GR_OFDM_ALAMOUTI_TX_CC_H
#define INCLUDED_GR_OFDM_ALAMOUTI_TX_CC_H

#include <gr_sync_block.h>
#include <vector>

class gr_ofdm_alamouti_tx_cc;
typedef boost::shared_ptr<gr_ofdm_alamouti_tx_cc> gr_ofdm_alamouti_tx_cc_sptr;

gr_ofdm_alamouti_tx_cc_sptr
gr_make_ofdm_alamouti_tx_cc(int fft_length);

/*!
 * \brief Take in OFDM symbols and perform Alamouti space-time coding for transmit
 * \ingroup ofdm_blk
 *
 * \param fft_length length of each symbol in samples.
 */

class gr_ofdm_alamouti_tx_cc : public gr_sync_block
{
  friend gr_ofdm_alamouti_tx_cc_sptr
  gr_make_ofdm_alamouti_tx_cc(int fft_length);

protected:
  gr_ofdm_alamouti_tx_cc(int fft_length);

private:
  int d_fft_length;

public:
  ~gr_ofdm_alamouti_tx_cc();

  int work (int noutput_items,
	    gr_vector_const_void_star &input_items,
	    gr_vector_void_star &output_items);
};

#endif /* INCLUDED_GR_OFDM_ALAMOUTI_TX_CC_H */
