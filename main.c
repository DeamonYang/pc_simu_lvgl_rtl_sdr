#define _DEFAULT_SOURCE /* needed for usleep() */
#include <stdlib.h>
#include <unistd.h>
#define SDL_MAIN_HANDLED /*To fix SDL's "undefined reference to WinMain" issue*/
#include <SDL2/SDL.h>
#include "lvgl/lvgl.h"
#include "lvgl/examples/lv_examples.h"
#include "lv_demos/lv_demo.h"
#include "lv_drivers/display/monitor.h"
#include "lv_drivers/display/fbdev.h"
#include "lv_drivers/indev/mouse.h"
#include "lv_drivers/indev/keyboard.h"
#include "lv_drivers/indev/mousewheel.h"
#include "rtl_sdr_app_work/include/rtl_sdr_dsp.h"
#include "tinyalsa/include/tinyalsa/asoundlib.h" 

#include<fftw3.h>
#include <errno.h>
#include <signal.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fftw3.h>
#include <math.h>
#include <pthread.h>
#include <libusb-1.0/libusb.h>

struct pcm *pcm;



#define PCM_PERIOD_SIZE 1024
#define PCM_PERIOD_COUNT 2

struct pcm_config  tiny_pcm_config = {
    .channels = 2,
    .format = PCM_FORMAT_S32_LE,
    .rate = 24000,
    .period_size = PCM_PERIOD_SIZE,
    .period_count = PCM_PERIOD_COUNT,
    .silence_threshold = PCM_PERIOD_SIZE*PCM_PERIOD_COUNT,
    .silence_size = 0,
    .stop_threshold = PCM_PERIOD_SIZE*PCM_PERIOD_COUNT,
    .start_threshold = PCM_PERIOD_SIZE
};

static volatile int do_exit = 0;
static int lcm_post[17] = {1,1,1,3,1,5,3,7,1,9,5,11,3,13,7,15,1};
static int ACTUAL_BUF_LENGTH;



static void hal_init(void);
static int tick_thread(void *data);
/*Pixel format: Fix 0xFF: 8 bit, Red: 8 bit, Green: 8 bit, Blue: 8 bit*/
//uint8_t disp_map[480*100*4];

lv_img_dsc_t fft_waterfall_dsc = {
  .header.always_zero = 0,
  .header.w = WATER_FALL_W,
  .header.h = WATER_FALL_H,
  .data_size = 48000 * LV_COLOR_SIZE / 8,
  .header.cf = LV_IMG_CF_TRUE_COLOR,
  //.data = disp_map,
};
lv_obj_t * img1 ;

static lv_obj_t * chart;
lv_chart_series_t * ser; 

// multiple of these, eventually
struct dongle_state dongle;
struct demod_state demod;
struct output_state output;
struct controller_state controller;
struct iq_fft_state iq_fft;



void pcm_player_init(struct pcm *pcm,struct pcm_config *config){
    pcm = pcm_open(0,0,PCM_OUT,config);
    if(!pcm_is_ready(pcm)){
        printf("failed to open pcm device\n");
        pcm_close(pcm);
    }
}

void pcm_play_load_buffer(struct pcm *pcm,void * data_buffer,int data_len){
    int write_frames = pcm_writei(pcm,data_buffer,pcm_bytes_to_frames(pcm,data_len));
    if(write_frames < 0){
        printf("error playing sample");
    } 
}


static void event_handler(lv_event_t * e)
{
    lv_event_code_t code = lv_event_get_code(e);
    lv_obj_t * obj = lv_event_get_target(e);
    if(code == LV_EVENT_VALUE_CHANGED) {
        char buf[32];
        lv_roller_get_selected_str(obj, buf, sizeof(buf));
        LV_LOG_USER("Selected value: %s", buf);
    }
}

void draw_init()
{
    img1= lv_img_create(lv_scr_act());
    lv_obj_set_pos (img1,  0, 172);
    lv_obj_set_size(img1,480, 100);    
    
    chart = lv_chart_create(lv_scr_act());
    lv_chart_set_type(chart, LV_CHART_TYPE_LINE);
    lv_obj_set_size(chart, 480, 100);
    lv_obj_set_pos(chart,0,72);    
        /*Do not display points on the data*/
    //lv_obj_set_style_size(chart, 0, LV_PART_INDICATOR);


    //lv_obj_add_event_cb(chart, draw_event_cb, LV_EVENT_DRAW_PART_BEGIN, NULL);
    //lv_chart_set_update_mode(chart, LV_CHART_UPDATE_MODE_CIRCULAR);
    lv_chart_set_point_count(chart, 480);
    lv_chart_set_div_line_count(chart, 10, 10);
    ser = lv_chart_add_series(chart, lv_palette_main(LV_PALETTE_LIGHT_BLUE), LV_CHART_AXIS_PRIMARY_Y);
    
    /*roller*/
    static lv_style_t style_sel;
    lv_style_init(&style_sel);
    lv_style_set_text_font(&style_sel, &lv_font_montserrat_22);

    const char * opts = "0\n1\n2\n3\n4\n5\n6\n7\n8\n9\n";
    const char * opts_po = ".0\n.1\n.2\n.3\n.4\n.5\n.6\n.7\n.8\n.9\n";
   
    lv_obj_t *rollerlist[10]; 
    for(int i = 0;i < 7;i++)
    {
        /*A roller on the left with left aligned text, and custom width*/
        rollerlist[i] = lv_roller_create(lv_scr_act());
        lv_obj_set_pos (rollerlist[i],  25*i, 0);
        lv_obj_set_size(rollerlist[i],25, 50);
        if(i == 4)
            lv_roller_set_options(rollerlist[i], opts_po, LV_ROLLER_MODE_NORMAL);
        else
            lv_roller_set_options(rollerlist[i], opts, LV_ROLLER_MODE_NORMAL);
        lv_roller_set_visible_row_count(rollerlist[i], 1);
        lv_obj_add_event_cb(rollerlist[i], event_handler, LV_EVENT_ALL, NULL);
        lv_roller_set_selected(rollerlist[i], i + 1 , LV_ANIM_OFF);  
    } 
 
    lv_obj_t * label2 = lv_label_create(lv_scr_act());
    lv_obj_set_pos (label2,  25*7, 10);
    lv_obj_set_size(label2,25*3, 50);
    //lv_label_set_long_mode(label2, LV_LABEL_LONG_SCROLL_CIRCULAR);     /*Circular scroll*/
    //lv_obj_set_width(label2, 150);
    lv_label_set_text(label2, "MHz");
    //lv_obj_align(label2, LV_ALIGN_CENTER, 0, 40); 
    
}

void draw_wave(lv_coord_t *ecg_sample,int max,int min)
{

    //lv_chart_set_range(chart, LV_CHART_AXIS_PRIMARY_Y, min, max);
    
    //lv_obj_add_event_cb(chart, draw_event_cb, LV_EVENT_DRAW_PART_BEGIN, NULL);
    //lv_chart_set_update_mode(chart, LV_CHART_UPDATE_MODE_CIRCULAR);
    lv_chart_set_ext_y_array(chart, ser, (lv_coord_t *)ecg_sample);        
}



static void sighandler(int signum)
{
	fprintf(stderr, "Signal caught, exiting!\n");
	do_exit = 1;
	rtlsdr_cancel_async(dongle.dev);
}


void fm_demod(struct demod_state *fm)
{
	int i, pcm;
	int16_t *lp = fm->lowpassed;
	pcm = polar_discriminant(lp[0], lp[1],
		fm->pre_r, fm->pre_j);
	fm->result[0] = (int16_t)pcm;
	for (i = 2; i < (fm->lp_len-1); i += 2) {
		switch (fm->custom_atan) {
		case 0:
			pcm = polar_discriminant(lp[i], lp[i+1],
				lp[i-2], lp[i-1]);
			break;
		case 1:
			pcm = polar_disc_fast(lp[i], lp[i+1],
				lp[i-2], lp[i-1]);
			break;
		case 2:
			pcm = polar_disc_lut(lp[i], lp[i+1],
				lp[i-2], lp[i-1]);
			break;
		}
		fm->result[i/2] = (int16_t)pcm;
	}
	fm->pre_r = lp[fm->lp_len - 2];
	fm->pre_j = lp[fm->lp_len - 1];
	fm->result_len = fm->lp_len/2;
}

void am_demod(struct demod_state *fm)
// todo, fix this extreme laziness
{
	int i, pcm;
	int16_t *lp = fm->lowpassed;
	int16_t *r  = fm->result;
	for (i = 0; i < fm->lp_len; i += 2) {
		// hypot uses floats but won't overflow
		//r[i/2] = (int16_t)hypot(lp[i], lp[i+1]);
		pcm = lp[i] * lp[i];
		pcm += lp[i+1] * lp[i+1];
		r[i/2] = (int16_t)sqrt(pcm) * fm->output_scale;
	}
	fm->result_len = fm->lp_len/2;
	// lowpass? (3khz)  highpass?  (dc)
}

void usb_demod(struct demod_state *fm)
{
	int i, pcm;
	int16_t *lp = fm->lowpassed;
	int16_t *r  = fm->result;
	for (i = 0; i < fm->lp_len; i += 2) {
		pcm = lp[i] + lp[i+1];
		r[i/2] = (int16_t)pcm * fm->output_scale;
	}
	fm->result_len = fm->lp_len/2;
}

void lsb_demod(struct demod_state *fm)
{
	int i, pcm;
	int16_t *lp = fm->lowpassed;
	int16_t *r  = fm->result;
	for (i = 0; i < fm->lp_len; i += 2) {
		pcm = lp[i] - lp[i+1];
		r[i/2] = (int16_t)pcm * fm->output_scale;
	}
	fm->result_len = fm->lp_len/2;
}

void raw_demod(struct demod_state *fm)
{
	int i;
	for (i = 0; i < fm->lp_len; i++) {
		fm->result[i] = (int16_t)fm->lowpassed[i];
	}
	fm->result_len = fm->lp_len;
}

void deemph_filter(struct demod_state *fm)
{
	static int avg;  // cheating...
	int i, d;
	// de-emph IIR
	// avg = avg * (1 - alpha) + sample * alpha;
	for (i = 0; i < fm->result_len; i++) {
		d = fm->result[i] - avg;
		if (d > 0) {
			avg += (d + fm->deemph_a/2) / fm->deemph_a;
		} else {
			avg += (d - fm->deemph_a/2) / fm->deemph_a;
		}
		fm->result[i] = (int16_t)avg;
	}
}

void dc_block_filter(struct demod_state *fm)
{
	int i, avg;
	int64_t sum = 0;
	for (i=0; i < fm->result_len; i++) {
		sum += fm->result[i];
	}
	avg = sum / fm->result_len;
	avg = (avg + fm->dc_avg * 9) / 10;
	for (i=0; i < fm->result_len; i++) {
		fm->result[i] -= avg;
	}
	fm->dc_avg = avg;
}


void full_demod(struct demod_state *d)
{
	int i, ds_p;
	int sr = 0;
	ds_p = d->downsample_passes;
	if (ds_p) {
		for (i=0; i < ds_p; i++) {
			fifth_order(d->lowpassed,   (d->lp_len >> i), d->lp_i_hist[i]);
			fifth_order(d->lowpassed+1, (d->lp_len >> i) - 1, d->lp_q_hist[i]);
		}
		d->lp_len = d->lp_len >> ds_p;
		/* droop compensation */
		if (d->comp_fir_size == 9 && ds_p <= CIC_TABLE_MAX) {
			generic_fir(d->lowpassed, d->lp_len,
				cic_9_tables[ds_p], d->droop_i_hist);
			generic_fir(d->lowpassed+1, d->lp_len-1,
				cic_9_tables[ds_p], d->droop_q_hist);
		}
	} else {
		low_pass(d);
	}
	/* power squelch */
	if (d->squelch_level) {
		sr = rms(d->lowpassed, d->lp_len, 1);
		if (sr < d->squelch_level) {
			d->squelch_hits++;
			for (i=0; i<d->lp_len; i++) {
				d->lowpassed[i] = 0;
			}
		} else {
			d->squelch_hits = 0;}
	}
	d->mode_demod(d);  /* lowpassed -> result */
	if (d->mode_demod == &raw_demod) {
		return;
	}
	/* todo, fm noise squelch */
	// use nicer filter here too?
	if (d->post_downsample > 1) {
		d->result_len = low_pass_simple(d->result, d->result_len, d->post_downsample);}
	if (d->deemph) {
		deemph_filter(d);}
	if (d->dc_block) {
		dc_block_filter(d);}
	if (d->rate_out2 > 0) {
		low_pass_real(d);
		//arbitrary_resample(d->result, d->result, d->result_len, d->result_len * d->rate_out2 / d->rate_out);
	}
}

static void rtlsdr_callback(unsigned char *buf, uint32_t len, void *ctx)
{
	int i;
	struct dongle_state *s = ctx;
	struct demod_state *d = s->demod_target;
	struct iq_fft_state *f = s->iq_fft_target;
	uint16_t buffer_temp[MAXIMUM_BUF_LENGTH];
	if (do_exit) {
		return;}
	if (!ctx) {
		return;}
	if (s->mute) {
		for (i=0; i<s->mute; i++) {
			buf[i] = 127;}
		s->mute = 0;
	}
	if (!s->offset_tuning) {
		rotate_90(buf, len);}
	for (i=0; i<(int)len; i++) {
		buffer_temp[i] = (int16_t)buf[i]-127;
        s->buf16[i] = (int16_t)buf[i]-127;//buffer_temp[i];
	}
	pthread_rwlock_wrlock(&d->rw);
	memcpy(d->lowpassed, s->buf16, 2*len);
	d->lp_len = len;
	pthread_rwlock_unlock(&d->rw);
	safe_cond_signal(&d->ready, &d->ready_m);
    /*copy data to fft buffer*/
    pthread_rwlock_wrlock(&f->rw);
    memcpy(f->c16buff,buffer_temp, 2*len);
    f->data_len = len;
	pthread_rwlock_unlock(&f->rw);
	safe_cond_signal(&f->ready, &f->ready_m);	
}

static void *dongle_thread_fn(void *arg)
{
	struct dongle_state *s = arg;
	rtlsdr_read_async(s->dev, rtlsdr_callback, s, 0, s->buf_len);
	return 0;
}

static void *demod_thread_fn(void *arg)
{
	struct demod_state *d = arg;
	struct output_state *o = d->output_target;
	while (!do_exit) {
		safe_cond_wait(&d->ready, &d->ready_m);
		pthread_rwlock_wrlock(&d->rw);
		full_demod(d);
		pthread_rwlock_unlock(&d->rw);
		if (d->exit_flag) {
			do_exit = 1;
		}
		if (d->squelch_level && d->squelch_hits > d->conseq_squelch) {
			d->squelch_hits = d->conseq_squelch + 1;  /* hair trigger */
			safe_cond_signal(&controller.hop, &controller.hop_m);
			continue;
		}
		pthread_rwlock_wrlock(&o->rw);
		memcpy(o->result, d->result, 2*d->result_len);
		o->result_len = d->result_len;
		pthread_rwlock_unlock(&o->rw);
		safe_cond_signal(&o->ready, &o->ready_m);
	}
	return 0;
}



static void *iq_fft_thread_fn(void *arg)
{
    struct iq_fft_state *f = arg;
    static int data_tot_len_acc = 0;
    int temp_data_len;
    int update_fft_num = 0;
    float fft_res_map[FFT_LEN];
    int i,j,k;
    int16_t data_buff[3][4096*2];
    fftw_complex in[FFT_LEN], out[FFT_LEN];
    fftw_plan p;
    float temp;
    int base_idx = 0;
    static unsigned int water_fall_map[WATER_FALL_W*WATER_FALL_H];
	static unsigned int (*ptmap)[WATER_FALL_W] = water_fall_map;
	float min,max;
	static lv_coord_t ecg_sample[FFT_LEN];
    
    while (!do_exit) {
        safe_cond_wait(&f->ready, &f->ready_m);
        pthread_rwlock_wrlock(&f->rw);
        update_fft_num = f->data_len/PR_FFT_LEN;
	update_fft_num = 1;
        if(update_fft_num > 0){
            for(i = 0; i< update_fft_num; i++){
                base_idx = i * PR_FFT_LEN;
                memcpy(data_buff[i],&f->c16buff[i * PR_FFT_LEN*2],FFT_LEN*sizeof(int16_t)*2);
            }
        }
        pthread_rwlock_unlock(&f->rw);


        for(i = 0; i< update_fft_num; i++){
            //memcpy(data_buff[i],&f->c16buff[i * PR_FFT_LEN*2],FFT_LEN*sizeof(int16_t)*2);
            for(j = 0;j < FFT_LEN;j++){
                in[j][0] = data_buff[i][j*2];
                in[j][1] = data_buff[i][j*2 + 1];
            }
	    
            
            p = fftw_plan_dft_1d(FFT_LEN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p); /* repeat as needed */
            
            for(j = 0;j < FFT_LEN;j ++){
                temp = out[j][0]*out[j][0] + out[j][1]*out[j][1];
                temp = sqrt(temp);
                fft_res_map[j] = 20*log10(temp);
            }
	    
		    // printf("fft_res_map[%d] data %0.2f %0.2f %0.2f\n",i,fft_res_map[0],fft_res_map[1],fft_res_map[2]);
            /******************waterfall & plot wave ***********************/
            max = 0;
			min = 9999999;
			for(j = 0;j < 480;j ++)
			{
				double temp = 0;
				for(k=0;k < 4;k++)
				{
					temp +=fft_res_map[j*8 + k];
				}
				
				ecg_sample[j] = temp/4;

				if(max < ecg_sample[j])
				{
					max = ecg_sample[j];
				}
				if(min > ecg_sample[j])
				{
					min = ecg_sample[j];
				}     
			}
#if 0
			printf("[max min][%d] is [%0.2f %0.2f]\n",i,max,min);
			FILE *fp;
			fp = fopen("dump_ecg_data.dat","w");
			for(j = 0;j < 480;j ++){
			    fprintf(fp,"%d ",ecg_sample[j]);
			}
			fprintf(fp,"\n");
			fclose(fp);
#endif
			draw_wave(ecg_sample,50,80);
				
			for(j = 0;j < WATER_FALL_H-1;j++){
				memcpy(ptmap[WATER_FALL_H-1-j],ptmap[WATER_FALL_H-2-j],WATER_FALL_W*sizeof(int));
			}
			
			for(j = 0;j < WATER_FALL_W;j++){
			    ptmap[0][j] = (unsigned int)ecg_sample[j]|0xFF000000|0x330000;
			}
			
			fft_waterfall_dsc.data=(uint8_t*)water_fall_map;
			lv_img_set_src(img1,&fft_waterfall_dsc);
			if(i < update_fft_num-1){
				;
				//usleep(10*1000);
			}
		}
    }
}




static void *output_thread_fn(void *arg)
{
	struct output_state *s = arg;
	while (!do_exit) {
		// use timedwait and pad out under runs
		safe_cond_wait(&s->ready, &s->ready_m);
		pthread_rwlock_rdlock(&s->rw);
		fwrite(s->result, 2, s->result_len, s->file);
		pthread_rwlock_unlock(&s->rw);
	}
	return 0;
}

static void optimal_settings(int freq, int rate)
{
	// giant ball of hacks
	// seems unable to do a single pass, 2:1
	int capture_freq, capture_rate;
	struct dongle_state *d = &dongle;
	struct demod_state *dm = &demod;
	struct controller_state *cs = &controller;
	dm->downsample = (1000000 / dm->rate_in) + 1;
	if (dm->downsample_passes) {
		dm->downsample_passes = (int)log2(dm->downsample) + 1;
		dm->downsample = 1 << dm->downsample_passes;
	}
	capture_freq = freq;
	capture_rate = dm->downsample * dm->rate_in;
	if (!d->offset_tuning) {
		capture_freq = freq + capture_rate/4;}
	capture_freq += cs->edge * dm->rate_in / 2;
	dm->output_scale = (1<<15) / (128 * dm->downsample);
	if (dm->output_scale < 1) {
		dm->output_scale = 1;}
	if (dm->mode_demod == &fm_demod) {
		dm->output_scale = 1;}
	d->freq = (uint32_t)capture_freq;
	d->rate = (uint32_t)capture_rate;
}

static void *controller_thread_fn(void *arg)
{
	// thoughts for multiple dongles
	// might be no good using a controller thread if retune/rate blocks
	int i;
	struct controller_state *s = arg;

	if (s->wb_mode) {
		for (i=0; i < s->freq_len; i++) {
			s->freqs[i] += 16000;}
	}

	/* set up primary channel */
	optimal_settings(s->freqs[0], demod.rate_in);
	if (dongle.direct_sampling) {
		verbose_direct_sampling(dongle.dev, 1);}
	if (dongle.offset_tuning) {
		verbose_offset_tuning(dongle.dev);}

	/* Set the frequency */
	verbose_set_frequency(dongle.dev, dongle.freq);
	fprintf(stderr, "Oversampling input by: %ix.\n", demod.downsample);
	fprintf(stderr, "Oversampling output by: %ix.\n", demod.post_downsample);
	fprintf(stderr, "Buffer size: %0.2fms\n",
		1000 * 0.5 * (float)ACTUAL_BUF_LENGTH / (float)dongle.rate);

	/* Set the sample rate */
	verbose_set_sample_rate(dongle.dev, dongle.rate);
	fprintf(stderr, "Output at %u Hz.\n", demod.rate_in/demod.post_downsample);

	while (!do_exit) {
		safe_cond_wait(&s->hop, &s->hop_m);
		if (s->freq_len <= 1) {
			continue;}
		/* hacky hopping */
		s->freq_now = (s->freq_now + 1) % s->freq_len;
		optimal_settings(s->freqs[s->freq_now], demod.rate_in);
		rtlsdr_set_center_freq(dongle.dev, dongle.freq);
		dongle.mute = BUFFER_DUMP;
	}
	return 0;
}

void frequency_range(struct controller_state *s, char *arg)
{
	char *start, *stop, *step;
	int i;
	start = arg;
	stop = strchr(start, ':') + 1;
	if (stop == (char *)1) { // no stop or step given
		s->freqs[s->freq_len] = (uint32_t) atofs(start);
		s->freq_len++;
		return;
	}
	stop[-1] = '\0';
	step = strchr(stop, ':') + 1;
	if (step == (char *)1) { // no step given
		s->freqs[s->freq_len] = (uint32_t) atofs(start);
		s->freq_len++;
		s->freqs[s->freq_len] = (uint32_t) atofs(stop);
		s->freq_len++;
		stop[-1] = ':';
		return;
	}
	step[-1] = '\0';
	for(i=(int)atofs(start); i<=(int)atofs(stop); i+=(int)atofs(step))
	{
		s->freqs[s->freq_len] = (uint32_t)i;
		s->freq_len++;
		if (s->freq_len >= FREQUENCIES_LIMIT) {
			break;}
	}
	stop[-1] = ':';
	step[-1] = ':';
}

void dongle_init(struct dongle_state *s)
{
	s->rate = DEFAULT_SAMPLE_RATE;
	s->gain = AUTO_GAIN; // tenths of a dB
	s->mute = 0;
	s->direct_sampling = 0;
	s->offset_tuning = 0;
	s->demod_target = &demod;
    s->iq_fft_target = &iq_fft;
}




void iq_fft_init(struct iq_fft_state *s)
{
    s->data_len = 0;
}

void demod_init(struct demod_state *s)
{
	s->rate_in = DEFAULT_SAMPLE_RATE;
	s->rate_out = DEFAULT_SAMPLE_RATE;
	s->squelch_level = 0;
	s->conseq_squelch = 10;
	s->terminate_on_squelch = 0;
	s->squelch_hits = 11;
	s->downsample_passes = 0;
	s->comp_fir_size = 0;
	s->prev_index = 0;
	s->post_downsample = 1;  // once this works, default = 4
	s->custom_atan = 0;
	s->deemph = 0;
	s->rate_out2 = -1;  // flag for disabled
	s->mode_demod = &fm_demod;
	s->pre_j = s->pre_r = s->now_r = s->now_j = 0;
	s->prev_lpr_index = 0;
	s->deemph_a = 0;
	s->now_lpr = 0;
	s->dc_block = 0;
	s->dc_avg = 0;
	pthread_rwlock_init(&s->rw, NULL);
	pthread_cond_init(&s->ready, NULL);
	pthread_mutex_init(&s->ready_m, NULL);
	s->output_target = &output;
}

void iq_fft_cleanup(struct iq_fft_state *s)
{
	pthread_rwlock_destroy(&s->rw);
	pthread_cond_destroy(&s->ready);
	pthread_mutex_destroy(&s->ready_m);
}

void demod_cleanup(struct demod_state *s)
{
	pthread_rwlock_destroy(&s->rw);
	pthread_cond_destroy(&s->ready);
	pthread_mutex_destroy(&s->ready_m);
}

void output_init(struct output_state *s)
{
	s->rate = DEFAULT_SAMPLE_RATE;
	pthread_rwlock_init(&s->rw, NULL);
	pthread_cond_init(&s->ready, NULL);
	pthread_mutex_init(&s->ready_m, NULL);
}

void output_cleanup(struct output_state *s)
{
	pthread_rwlock_destroy(&s->rw);
	pthread_cond_destroy(&s->ready);
	pthread_mutex_destroy(&s->ready_m);
}

void controller_init(struct controller_state *s)
{
	s->freqs[0] = 100000000;
	s->freq_len = 0;
	s->edge = 0;
	s->wb_mode = 0;
	pthread_cond_init(&s->hop, NULL);
	pthread_mutex_init(&s->hop_m, NULL);
}

void controller_cleanup(struct controller_state *s)
{
	pthread_cond_destroy(&s->hop);
	pthread_mutex_destroy(&s->hop_m);
}

void sanity_checks(void)
{
	if (controller.freq_len == 0) {
		fprintf(stderr, "Please specify a frequency.\n");
		exit(1);
	}

	if (controller.freq_len >= FREQUENCIES_LIMIT) {
		fprintf(stderr, "Too many channels, maximum %i.\n", FREQUENCIES_LIMIT);
		exit(1);
	}

	if (controller.freq_len > 1 && demod.squelch_level == 0) {
		fprintf(stderr, "Please specify a squelch level.  Required for scanning multiple frequencies.\n");
		exit(1);
	}

}



int main(int argc, char **argv)
{
	struct sigaction sigact;
	int r, opt;
	int dev_given = 0;
	int custom_ppm = 0;
    int enable_biastee = 0;
	lv_timer_t *timer_waterfull;
	dongle_init(&dongle);
	demod_init(&demod);
	output_init(&output);
	controller_init(&controller);
    iq_fft_init(&iq_fft);
    pcm_player_init(pcm,&tiny_pcm_config);
	
	
    lv_init();
    hal_init();
    draw_init();
    //timer_waterfull = lv_timer_create(draw_fft_waterfull,100,0); //period_ms, user_data

	while ((opt = getopt(argc, argv, "d:f:g:s:b:l:o:t:r:p:E:F:A:M:hT")) != -1) {
		switch (opt) {
		case 'd':
			dongle.dev_index = verbose_device_search(optarg);
			dev_given = 1;
			break;
		case 'f':
			if (controller.freq_len >= FREQUENCIES_LIMIT) {
				break;}
			if (strchr(optarg, ':'))
				{frequency_range(&controller, optarg);}
			else
			{
				controller.freqs[controller.freq_len] = (uint32_t)atofs(optarg);
				controller.freq_len++;
			}
			break;
		case 'g':
			dongle.gain = (int)(atof(optarg) * 10);
			break;
		case 'l':
			demod.squelch_level = (int)atof(optarg);
			break;
		case 's':
			demod.rate_in = (uint32_t)atofs(optarg);
			demod.rate_out = (uint32_t)atofs(optarg);
			break;
		case 'r':
			output.rate = (int)atofs(optarg);
			demod.rate_out2 = (int)atofs(optarg);
			break;
		case 'o':
			fprintf(stderr, "Warning: -o is very buggy\n");
			demod.post_downsample = (int)atof(optarg);
			if (demod.post_downsample < 1 || demod.post_downsample > MAXIMUM_OVERSAMPLE) {
				fprintf(stderr, "Oversample must be between 1 and %i\n", MAXIMUM_OVERSAMPLE);}
			break;
		case 't':
			demod.conseq_squelch = (int)atof(optarg);
			if (demod.conseq_squelch < 0) {
				demod.conseq_squelch = -demod.conseq_squelch;
				demod.terminate_on_squelch = 1;
			}
			break;
		case 'p':
			dongle.ppm_error = atoi(optarg);
			custom_ppm = 1;
			break;
		case 'E':
			if (strcmp("edge",  optarg) == 0) {
				controller.edge = 1;}
			if (strcmp("dc", optarg) == 0) {
				demod.dc_block = 1;}
			if (strcmp("deemp",  optarg) == 0) {
				demod.deemph = 1;}
			if (strcmp("direct",  optarg) == 0) {
				dongle.direct_sampling = 1;}
			if (strcmp("offset",  optarg) == 0) {
				dongle.offset_tuning = 1;}
			break;
		case 'F':
			demod.downsample_passes = 1;  /* truthy placeholder */
			demod.comp_fir_size = atoi(optarg);
			break;
		case 'A':
			if (strcmp("std",  optarg) == 0) {
				demod.custom_atan = 0;}
			if (strcmp("fast", optarg) == 0) {
				demod.custom_atan = 1;}
			if (strcmp("lut",  optarg) == 0) {
				atan_lut_init();
				demod.custom_atan = 2;}
			break;
		case 'M':
			if (strcmp("fm",  optarg) == 0) {
				demod.mode_demod = &fm_demod;}
			if (strcmp("raw",  optarg) == 0) {
				demod.mode_demod = &raw_demod;}
			if (strcmp("am",  optarg) == 0) {
				demod.mode_demod = &am_demod;}
			if (strcmp("usb", optarg) == 0) {
				demod.mode_demod = &usb_demod;}
			if (strcmp("lsb", optarg) == 0) {
				demod.mode_demod = &lsb_demod;}
			if (strcmp("wbfm",  optarg) == 0) {
				controller.wb_mode = 1;
				demod.mode_demod = &fm_demod;
				demod.rate_in = 170000;
				demod.rate_out = 170000;
				demod.rate_out2 = 32000;
				demod.custom_atan = 1;
				//demod.post_downsample = 4;
				demod.deemph = 1;
				demod.squelch_level = 0;}
			break;
		case 'T':
			enable_biastee = 1;
			break;
		case 'h':
		default:
			usage(); 
			break;
		}
	}

	/* quadruple sample_rate to limit to Δθ to ±π/2 */
	demod.rate_in *= demod.post_downsample;

	if (!output.rate) {
		output.rate = demod.rate_out;}

	sanity_checks();

	if (controller.freq_len > 1) {
		demod.terminate_on_squelch = 0;}

	if (argc <= optind) {
		output.filename = "-";
	} else {
		output.filename = argv[optind];
	}

	ACTUAL_BUF_LENGTH = lcm_post[demod.post_downsample] * DEFAULT_BUF_LENGTH;

	if (!dev_given) {
		dongle.dev_index = verbose_device_search("0");
	}

	if (dongle.dev_index < 0) {
		exit(1);
	}

	r = rtlsdr_open(&dongle.dev, (uint32_t)dongle.dev_index);
	if (r < 0) {
		fprintf(stderr, "Failed to open rtlsdr device #%d.\n", dongle.dev_index);
		exit(1);
	}

	sigact.sa_handler = sighandler;
	sigemptyset(&sigact.sa_mask);
	sigact.sa_flags = 0;
	sigaction(SIGINT, &sigact, NULL);
	sigaction(SIGTERM, &sigact, NULL);
	sigaction(SIGQUIT, &sigact, NULL);
	sigaction(SIGPIPE, &sigact, NULL);


	if (demod.deemph) {
		demod.deemph_a = (int)round(1.0/((1.0-exp(-1.0/(demod.rate_out * 75e-6)))));
	}

	/* Set the tuner gain */
	if (dongle.gain == AUTO_GAIN) {
		verbose_auto_gain(dongle.dev);
	} else {
		dongle.gain = nearest_gain(dongle.dev, dongle.gain);
		verbose_gain_set(dongle.dev, dongle.gain);
	}

	rtlsdr_set_bias_tee(dongle.dev, enable_biastee);
	if (enable_biastee)
		fprintf(stderr, "activated bias-T on GPIO PIN 0\n");

	verbose_ppm_set(dongle.dev, dongle.ppm_error);

	if (strcmp(output.filename, "-") == 0) { /* Write samples to stdout */
		output.file = stdout;
	} else {
		output.file = fopen(output.filename, "wb");
		if (!output.file) {
			fprintf(stderr, "Failed to open %s\n", output.filename);
			exit(1);
		}
	}

	//r = rtlsdr_set_testmode(dongle.dev, 1);

	/* Reset endpoint before we start reading from it (mandatory) */
	verbose_reset_buffer(dongle.dev);

	pthread_create(&controller.thread, NULL, controller_thread_fn, (void *)(&controller));
	usleep(100000);
	pthread_create(&output.thread, NULL, output_thread_fn, (void *)(&output));
    pthread_create(&iq_fft.thread, NULL, iq_fft_thread_fn, (void *)(&iq_fft));
	pthread_create(&demod.thread, NULL, demod_thread_fn, (void *)(&demod));
	pthread_create(&dongle.thread, NULL, dongle_thread_fn, (void *)(&dongle));
	
    while (!do_exit) {
	    lv_timer_handler();
		usleep(5 * 1000);
	}	
	
	if (do_exit) {
		fprintf(stderr, "\nUser cancel, exiting...\n");}
	else {
		fprintf(stderr, "\nLibrary error %d, exiting...\n", r);}

	rtlsdr_cancel_async(dongle.dev);
	pthread_join(dongle.thread, NULL);
    
	safe_cond_signal(&demod.ready, &demod.ready_m);
	pthread_join(demod.thread, NULL);
    
    safe_cond_signal(&iq_fft.ready, &iq_fft.ready_m);
    pthread_join(iq_fft.thread, NULL);
    
	safe_cond_signal(&output.ready, &output.ready_m);
	pthread_join(output.thread, NULL);
    
	safe_cond_signal(&controller.hop, &controller.hop_m);
	pthread_join(controller.thread, NULL);
   

	//dongle_cleanup(&dongle);
	demod_cleanup(&demod);
    iq_fft_cleanup(&iq_fft);
	output_cleanup(&output);
	controller_cleanup(&controller);

	if (output.file != stdout) {
		fclose(output.file);}

	rtlsdr_close(dongle.dev);
	return r >= 0 ? r : -r;

	
}



/**********************
 *   STATIC FUNCTIONS
 **********************/

/**
 * Initialize the Hardware Abstraction Layer (HAL) for the LVGL graphics
 * library
 */
static void hal_init(void)
{
  /* Use the 'monitor' driver which creates window on PC's monitor to simulate a display*/
  monitor_init();
  /* Tick init.
   * You have to call 'lv_tick_inc()' in periodically to inform LittelvGL about
   * how much time were elapsed Create an SDL thread to do this*/
  SDL_CreateThread(tick_thread, "tick", NULL);

  /*Create a display buffer*/
  static lv_disp_draw_buf_t disp_buf1;
  static lv_color_t buf1_1[MONITOR_HOR_RES * 100];
  static lv_color_t buf1_2[MONITOR_HOR_RES * 100];
  lv_disp_draw_buf_init(&disp_buf1, buf1_1, buf1_2, MONITOR_HOR_RES * 100);

  /*Create a display*/
  static lv_disp_drv_t disp_drv;
  lv_disp_drv_init(&disp_drv); /*Basic initialization*/
  disp_drv.draw_buf = &disp_buf1;
  disp_drv.flush_cb = monitor_flush;
  disp_drv.hor_res = MONITOR_HOR_RES;
  disp_drv.ver_res = MONITOR_VER_RES;
  disp_drv.antialiasing = 1;

  lv_disp_t * disp = lv_disp_drv_register(&disp_drv);

  lv_theme_t * th = lv_theme_default_init(disp, lv_palette_main(LV_PALETTE_BLUE), lv_palette_main(LV_PALETTE_RED), LV_THEME_DEFAULT_DARK, LV_FONT_DEFAULT);
  lv_disp_set_theme(disp, th);

  lv_group_t * g = lv_group_create();
  lv_group_set_default(g);

  /* Add the mouse as input device
   * Use the 'mouse' driver which reads the PC's mouse*/
  mouse_init();
  static lv_indev_drv_t indev_drv_1;
  lv_indev_drv_init(&indev_drv_1); /*Basic initialization*/
  indev_drv_1.type = LV_INDEV_TYPE_POINTER;

  /*This function will be called periodically (by the library) to get the mouse position and state*/
  indev_drv_1.read_cb = mouse_read;
  lv_indev_t *mouse_indev = lv_indev_drv_register(&indev_drv_1);

  keyboard_init();
  static lv_indev_drv_t indev_drv_2;
  lv_indev_drv_init(&indev_drv_2); /*Basic initialization*/
  indev_drv_2.type = LV_INDEV_TYPE_KEYPAD;
  indev_drv_2.read_cb = keyboard_read;
  lv_indev_t *kb_indev = lv_indev_drv_register(&indev_drv_2);
  lv_indev_set_group(kb_indev, g);
  mousewheel_init();
  static lv_indev_drv_t indev_drv_3;
  lv_indev_drv_init(&indev_drv_3); /*Basic initialization*/
  indev_drv_3.type = LV_INDEV_TYPE_ENCODER;
  indev_drv_3.read_cb = mousewheel_read;

  lv_indev_t * enc_indev = lv_indev_drv_register(&indev_drv_3);
  lv_indev_set_group(enc_indev, g);

  /*Set a cursor for the mouse*/
  LV_IMG_DECLARE(mouse_cursor_icon); /*Declare the image file.*/
  lv_obj_t * cursor_obj = lv_img_create(lv_scr_act()); /*Create an image object for the cursor */
  lv_img_set_src(cursor_obj, &mouse_cursor_icon);           /*Set the image source*/
  lv_indev_set_cursor(mouse_indev, cursor_obj);             /*Connect the image  object to the driver*/
}

/**
 * A task to measure the elapsed time for LVGL
 * @param data unused
 * @return never return
 */
static int tick_thread(void *data) {
    (void)data;

    while(1) { 
        SDL_Delay(30);
        lv_tick_inc(30); /*Tell LittelvGL that 5 milliseconds were elapsed*/
    }

    return 0;
}


// 
