#!/usr/bin/env python
#
# Copyright 2006, 2007, 2008 Free Software Foundation, Inc.
# 
# This file is part of GNU Radio
# 
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

import math
from numpy import fft
from gnuradio import gr
from gnuradio.blks2impl.ofdm_sync_ml import ofdm_sync_ml
from gnuradio.blks2impl.ofdm_sync_pn import ofdm_sync_pn
from gnuradio.blks2impl.ofdm_sync_pnac import ofdm_sync_pnac
from gnuradio.blks2impl.ofdm_sync_fixed import ofdm_sync_fixed

class ofdm_mimo_receiver(gr.hier_block2):
    """
    Performs receiver synchronization on OFDM symbols.

    The receiver performs channel filtering as well as symbol, frequency, and phase synchronization.
    The synchronization routines are available in three flavors: preamble correlator (Schmidl and Cox),
    modifid preamble correlator with autocorrelation (not yet working), and cyclic prefix correlator
    (Van de Beeks).
    """

    def __init__(self, ntxchans, nrxchans, fft_length, cp_length, occupied_tones, snr, ks, logging=False):
        """
	Hierarchical block for receiving OFDM symbols.

	The input is the complex modulated signal at baseband.
        Synchronized packets are sent back to the demodulator.

        @param ntxchans: Number of transmitter MIMO channels (antennas)
        @type  ntxchans: int
        @param nrxchans: Number of reciever MIMO channels (antennas)
        @type  nrxchans: int
        @param fft_length: total number of subcarriers
        @type  fft_length: int
        @param cp_length: length of cyclic prefix as specified in subcarriers (<= fft_length)
        @type  cp_length: int
        @param occupied_tones: number of subcarriers used for data
        @type  occupied_tones: int
        @param snr: estimated signal to noise ratio used to guide cyclic prefix synchronizer
        @type  snr: float
        @param ks: known symbols used as preambles to each packet
        @type  ks: list of lists
        @param logging: turn file logging on or off
        @type  logging: bool
	"""

	gr.hier_block2.__init__(self, "ofdm_mimo_receiver",
				gr.io_signature(1, 1, gr.sizeof_gr_complex), # Input signature
                                gr.io_signature2(2, 2, gr.sizeof_gr_complex*occupied_tones, gr.sizeof_char)) # Output signature
        
        bw = (float(occupied_tones) / float(fft_length)) / 2.0
        tb = bw*0.08
        chan_coeffs = gr.firdes.low_pass (1.0,                     # gain
                                          1.0,                     # sampling rate
                                          bw+tb,                   # midpoint of trans. band
                                          tb,                      # width of trans. band
                                          gr.firdes.WIN_HAMMING)   # filter type

        # For starters, run Sync on channel 0 and use it to clock and retime all channels
        win = [1 for i in range(fft_length)]

        zeros_on_left = int(math.ceil((fft_length - occupied_tones)/2.0))
        ks0 = fft_length*[0,]
        ks0[zeros_on_left : zeros_on_left + occupied_tones] = ks[0][0:]
        
        ks0 = fft.ifftshift(ks0)
        ks0time = fft.ifft(ks0)
        # ADD SCALING FACTOR
        ks0time = ks0time.tolist()

        SYNC = "pn"
        if SYNC == "ml":
            nco_sensitivity = -1.0/fft_length                             # correct for fine frequency
            self.ofdm_sync = ofdm_sync_ml(fft_length, cp_length, snr, ks0time, logging)
        elif SYNC == "pn":
            nco_sensitivity = -2.0/fft_length                             # correct for fine frequency
            self.ofdm_sync = ofdm_sync_pn(fft_length, cp_length, logging)
        elif SYNC == "pnac":
            nco_sensitivity = -2.0/fft_length                             # correct for fine frequency
            self.ofdm_sync = ofdm_sync_pnac(fft_length, cp_length, ks0time, logging)
        elif SYNC == "fixed":                                             # for testing only; do not user over the air
            self.chan_filt = gr.multiply_const_cc(1.0)                    # remove filter and filter delay for this
            nsymbols = 18                                                 # enter the number of symbols per packet
            freq_offset = 0.0                                             # if you use a frequency offset, enter it here
            nco_sensitivity = -2.0/fft_length                             # correct for fine frequency
            self.ofdm_sync = ofdm_sync_fixed(fft_length, cp_length, nsymbols, freq_offset, logging)

        # Set up blocks
        self.chan_filt = list()
        self.sigmix = list()
        self.sampler = list()
        self.fft_demod = list()

        # Deinterleave the incoming stream into separate channels
        self.deint = gr.deinterleave(gr.sizeof_gr_complex)

        # generate a signal proportional to frequency error of sync block
        self.nco = gr.frequency_modulator_fc(nco_sensitivity)

        # Manage and combine all channels
        if 0:
            self.ofdm_frame_acq = gr.ofdm_mrc_frame_acquisition(nrxchans, occupied_tones, fft_length,
                                                                cp_length, ks[0])
        else:
            print "Transmitter using ", ntxchans
            print "Receiver using    ", nrxchans
            if(ntxchans == 1):
                self.ofdm_frame_acq = gr.ofdm_alamouti_frame_acquisition(ntxchans, occupied_tones, fft_length,
                                                                         cp_length, ks[0], occupied_tones*[0,])
            else:
                self.ofdm_frame_acq = gr.ofdm_alamouti_frame_acquisition(ntxchans, occupied_tones, fft_length,
                                                                         cp_length, ks[0], ks[1])

        self.connect(self, self.deint)               # deinterleave channels
        self.connect((self.ofdm_sync,0), self.nco)   # use sync freq. offset to derotate signal

        for i in xrange(nrxchans):
            self.chan_filt.append(gr.fft_filter_ccc(1, chan_coeffs))
            self.sigmix.append(gr.multiply_cc())
            self.sampler.append(gr.ofdm_sampler(fft_length, fft_length+cp_length))
            self.fft_demod.append(gr.fft_vcc(fft_length, True, win, True))
            

            self.connect((self.deint, i), self.chan_filt[i])              # filter the input channel
            self.connect(self.nco, (self.sigmix[i],1))                    # use sync freq. offset to derotate signal
            self.connect(self.chan_filt[i], (self.sigmix[i],0))           # signal to be derotated
            self.connect(self.sigmix[i], (self.sampler[i],0))             # sample off timing signal detected in sync alg
            self.connect((self.ofdm_sync,1), (self.sampler[i],1))         # timing signal to sample at
            
            self.connect((self.sampler[i],0), self.fft_demod[i])          # send derotated sampled signal to FFT
            self.connect(self.fft_demod[i], (self.ofdm_frame_acq,1+i))    # find frame start and equalize signal

            if logging:
                self.connect(self.chan_filt[i],
                             gr.file_sink(gr.sizeof_gr_complex, ("ofdm_mimo-receiver-chan%02d-chan_filt_c.dat" % i)))
                self.connect(self.fft_demod[i],
                             gr.file_sink(gr.sizeof_gr_complex*fft_length, ("ofdm_mimo-receiver-chan%02d-fft_out_c.dat" % i)))
                self.connect(self.sampler[i],
                             gr.file_sink(gr.sizeof_gr_complex*fft_length, ("ofdm_mimo-receiver-chan%02d-sampler_c.dat" % i)))
                self.connect(self.sigmix[i],
                             gr.file_sink(gr.sizeof_gr_complex, ("ofdm_mimo-receiver-chan%02d-sigmix_c.dat" % i)))
                
        if logging:
            self.connect((self.ofdm_frame_acq,0),
                         gr.file_sink(gr.sizeof_gr_complex*occupied_tones, "ofdm_mimo-receiver-frame_acq_c.dat"))
            self.connect((self.ofdm_frame_acq,1), gr.file_sink(1, "ofdm_mimo-receiver-found_corr_b.dat"))
            self.connect(self.nco, gr.file_sink(gr.sizeof_gr_complex, "ofdm_mimo-receiver-nco_c.dat"))
            
        self.connect(self.chan_filt[0], self.ofdm_sync)               # into the synchronization alg.
        self.connect((self.sampler[0],1), (self.ofdm_frame_acq,0))    # send timing signal to signal frame start
        self.connect((self.sampler[i],1), gr.null_sink(fft_length*gr.sizeof_char))
        self.connect((self.ofdm_frame_acq,0), (self,0))               # finished with fine/coarse freq correction,
        self.connect((self.ofdm_frame_acq,1), (self,1))               # frame and symbol timing, and equalization
