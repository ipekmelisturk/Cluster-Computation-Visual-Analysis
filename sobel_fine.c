/* Example sobel code for ECE574 -- Sample HW6 Spring 2023 */
/* By Vince Weaver <vincent.weaver@maine.edu> */

/* NOTE: this code does not properly handle cases where ysize */
/*       is not a multiple of rank number */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <math.h>

#include <jpeglib.h>

#include <mpi.h>

/* Filters */
static int sobel_x_filter[9]={ -1, 0,+1  ,-2, 0,+2  ,-1, 0,+1};
static int sobel_y_filter[9]={ -1,-2,-1  , 0, 0, 0  , 1, 2,+1};

/* Structure describing the image */
struct image_t {
	int xsize;
	int ysize;
	int depth;	/* bytes */
	unsigned char *pixels;
};

struct convolve_data_t {
        struct image_t *old;
        struct image_t *new;
        int *filter;
        int ystart;
        int yend;
};


/* very inefficient convolve code */
static void generic_convolve(struct convolve_data_t *data) {

	int x,y,c,k,l;
	int color;
	int sum;

	int y_start = 0;
	int y_end = 0;

	struct image_t *input;
        struct image_t *output;
	int *filter;

	input=data->old;
	output=data->new;

	filter=data->filter;

	/* Don't run sobel on first or last row */
	if (data->ystart == 0) {
		y_start = 1;
	}
	else {
		y_start=data->ystart;
	}

	if (data->yend >= input->ysize) {
		y_end = input->ysize - 1;
	}
	else {
		y_end = data->yend;
	}

//	printf("Sobel from %d to %d\n",data->ystart,data->yend);
//	printf("\tAdjusted: %d to %d\n",y_start,y_end);

	for(c=0;c<(input->depth);c++) {
	   for(x=1;x<(input->xsize)-1;x++) {
              for(y=y_start;y<y_end;y++) {

		sum=0;

		for(k=-1;k<2;k++) {
		   for(l=-1;l<2;l++) {
			color=input->pixels[((y+l)*input->xsize*input->depth)+
                                            ((x+k)*input->depth)+
						c];
                        sum+=color * filter[(l+1)*3+(k+1)];
                   }
                }

		/* saturate */
		if (sum<0) sum=0;
		if (sum>255) sum=255;

		/* set output pixel */
		/* NOTE: we subtract off row_start so result */
		/*       is at the beginning of array */
		/*       this is needed when using plain MPI_Gather() */
		output->pixels[(y - data->ystart)*input->xsize*input->depth+
				x*input->depth
				+c]=sum;
	     }
	   }
	}
}


static int combine(struct image_t *s_x,
			struct image_t *s_y,
			struct image_t *new) {
	int i;
	int out;

	for(i=0;i<( s_x->depth * s_x->xsize * s_x->ysize );i++) {

		out=sqrt(
			(s_x->pixels[i]*s_x->pixels[i])+
			(s_y->pixels[i]*s_y->pixels[i])
			);
		if (out>255) out=255;
		if (out<0) out=0;
		new->pixels[i]=out;
	}

	return 0;
}



static int load_jpeg(char *filename, struct image_t *image) {

	FILE *fff;
	struct jpeg_decompress_struct cinfo;
	struct jpeg_error_mgr jerr;
	JSAMPROW output_data;
	unsigned int scanline_len;
	int scanline_count=0;

	fff=fopen(filename,"rb");
	if (fff==NULL) {
		fprintf(stderr, "Could not load %s: %s\n",
			filename, strerror(errno));
		return -1;
	}

	/* set up jpeg error routines */
	cinfo.err = jpeg_std_error(&jerr);

	/* Initialize cinfo */
	jpeg_create_decompress(&cinfo);

	/* Set input file */
	jpeg_stdio_src(&cinfo, fff);

	/* read header */
	jpeg_read_header(&cinfo, TRUE);

	/* Start decompressor */
	jpeg_start_decompress(&cinfo);

	printf("output_width=%d, output_height=%d, output_components=%d\n",
		cinfo.output_width,
		cinfo.output_height,
		cinfo.output_components);

	image->xsize=cinfo.output_width;
	image->ysize=cinfo.output_height;
	image->depth=cinfo.output_components;

	scanline_len = cinfo.output_width * cinfo.output_components;
	image->pixels=malloc(cinfo.output_width * cinfo.output_height * cinfo.output_components);

	while (scanline_count < cinfo.output_height) {
		output_data = (image->pixels + (scanline_count * scanline_len));
		jpeg_read_scanlines(&cinfo, &output_data, 1);
		scanline_count++;
	}

	/* Finish decompressing */
	jpeg_finish_decompress(&cinfo);

	jpeg_destroy_decompress(&cinfo);

	fclose(fff);

	return 0;
}

static int store_jpeg(char *filename, struct image_t *image) {

	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr jerr;
	int quality=90; /* % */
	int i;

	FILE *fff;

	JSAMPROW row_pointer[1];
	int row_stride;

	/* setup error handler */
	cinfo.err = jpeg_std_error(&jerr);

	/* initialize jpeg compression object */
	jpeg_create_compress(&cinfo);

	/* Open file */
	fff = fopen(filename, "wb");
	if (fff==NULL) {
		fprintf(stderr, "can't open %s: %s\n",
			filename,strerror(errno));
		return -1;
	}

	jpeg_stdio_dest(&cinfo, fff);

	/* Set compression parameters */
	cinfo.image_width = image->xsize;
	cinfo.image_height = image->ysize;
	cinfo.input_components = image->depth;
	cinfo.in_color_space = JCS_RGB;
	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, quality, TRUE);

	/* start compressing */
	jpeg_start_compress(&cinfo, TRUE);

	row_stride=image->xsize*image->depth;

	for(i=0;i<image->ysize;i++) {
		row_pointer[0] = & image->pixels[i * row_stride];
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	/* finish compressing */
	jpeg_finish_compress(&cinfo);

	/* close file */
	fclose(fff);

	/* clean up */
	jpeg_destroy_compress(&cinfo);

	return 0;
}


int main(int argc, char **argv) {

	struct image_t image,sobel_x,sobel_y,new_image;
	struct image_t result_x,result_y;
	double start_time,load_time,store_time,convolve_time,combine_time;
	int result, image_dim[3];
    int rank, size_mem, num_ranks;
	int row_start,row_end;
	int gather_size;
	int *displs,*recv_counts;
	//int *image_sc;
	struct convolve_data_t convolve_data;

//	MPI_Status stat;

	/* Check command line usage */
	if (argc<2) {
		fprintf(stderr,"Usage: %s image_file\n",argv[0]);
		return -1;
	}
        /* Initialize MPI */
	MPI_Init(&argc,&argv);

	// call rank and number of processors
	start_time=MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


        /* Load the image, only in Rank0 */
	if (rank==0) {
		result=load_jpeg(argv[1],&image);

		if (result<0) {
			fprintf(stderr,"Unable to load image %s\n",argv[1]);
			return -1;
		}


		/* get image data ready to broadcast */
		image_dim[0]= image.xsize;
		image_dim[1]= image.ysize;
		image_dim[2] = image.depth;
	}

	load_time=MPI_Wtime();

	/* Broadcast image sizes so that we can use them on all ranks */
	result = MPI_Bcast(image_dim,
				3,
				MPI_INT,
				0,
				MPI_COMM_WORLD);

	/* Now image.dim[] should be same on all ranks */

	size_mem = image_dim[0]*image_dim[1]*image_dim[2]*sizeof(char);

	sobel_x.xsize=image_dim[0];
	sobel_x.ysize=image_dim[1];
	sobel_x.depth=image_dim[2];
	sobel_x.pixels=calloc(size_mem,sizeof(char));

	sobel_y.xsize=image_dim[0];
	sobel_y.ysize=image_dim[1];
	sobel_y.depth=image_dim[2];
	sobel_y.pixels=calloc(size_mem,sizeof(char));

	result_x.xsize=image_dim[0];
	result_x.ysize=image_dim[1];
	result_x.depth=image_dim[2];
	result_x.pixels=calloc(size_mem,sizeof(char));

	result_y.xsize=image_dim[0];
	result_y.ysize=image_dim[1];
	result_y.depth=image_dim[2];
	result_y.pixels=calloc(size_mem,sizeof(char));

	new_image.xsize=image_dim[0];
	new_image.ysize=image_dim[1];
	new_image.depth=image_dim[2];
	new_image.pixels=calloc(size_mem,sizeof(char));

	/* Allocate image structs on other ranks */
	if (rank != 0) {
		image.xsize = sobel_x.xsize;
		image.ysize = sobel_x.ysize;
        image.depth = sobel_x.depth;
		image.pixels = calloc(size_mem,sizeof(char));
	}

	
	result = MPI_Bcast(image.pixels,
		image.xsize*image.ysize*image.depth,
		MPI_CHAR,
		0,
		MPI_COMM_WORLD);


	/* Setup row split */

	row_start=rank*(sobel_x.ysize/num_ranks);

	// Set the row_end to the reaminder if it is the last rank 
	if (rank == num_ranks)
	{
		row_end = sobel_x.ysize;
	}
	else
	{
		row_end=(rank+1)*(sobel_x.ysize/num_ranks);
	}
	
	

	//printf("Rank %d: calculating from %d to %d\n",rank,row_start,row_end);

	/***********************/
	/* Sobel X convolution */
	/***********************/

	convolve_data.old=&image;
	convolve_data.new=&sobel_x;
	convolve_data.filter=sobel_x_filter;

	/* Setup row split */
    convolve_data.ystart=row_start;
	convolve_data.yend=row_end;
	printf("SobelX Rank %d: calculating from %d to %d\n",rank,
	convolve_data.ystart,
	convolve_data.yend);

	generic_convolve(&convolve_data);
	/***********************/
	/* Gather in Y Results */ 
	/***********************/
	
	// For Gatherv results
	
	displs = (int *)calloc(size_mem,sizeof(char));
    recv_counts = (int *)calloc(size_mem,sizeof(char));

	for(int i=0; i<num_ranks; ++i )
	{
		recv_counts[i] = (int *)((sobel_x.ysize/num_ranks) * sobel_x.xsize*sobel_x.depth) ;
		displs[i] = i*(sobel_x.ysize/num_ranks) * sobel_x.xsize*sobel_x.depth  ;
			
	}
	

	gather_size=(sobel_x.ysize/num_ranks)*sobel_x.xsize*sobel_x.depth;
	if (rank==0) {
		printf("Using gather size of %d\n",gather_size);
	}

	/* Gather the data from the start of each rank's sobel_x.pixels */
	/* and put it at the proper offset in rank 0's result_x.pixels */

	/* Note there's an implicit barrier here */
	
		MPI_Gatherv(sobel_x.pixels,
		gather_size,
		MPI_CHAR,
		result_x.pixels,
		recv_counts,
		displs,
		MPI_CHAR,
		0,
		MPI_COMM_WORLD);


	/***********************/
	/* Sobel Y convolution */
	/***********************/

	convolve_data.old=&image;
	convolve_data.new=&sobel_y;
	convolve_data.filter=sobel_y_filter;

	/* Setup row split */
      convolve_data.ystart=rank*(sobel_y.ysize/num_ranks);
	convolve_data.yend=(rank+1)*(sobel_y.ysize/num_ranks);
	printf("SobelY Rank %d: calculating from %d to %d\n",rank,
		convolve_data.ystart,
		convolve_data.yend);

	/* Do the Y convolution */
	generic_convolve(&convolve_data);

	/***********************/
	/* Gather in X Results */
	/***********************/

	/* Gather the data from the start of each rank's sobel_y.pixels */
	/* and put it at the proper offset in rank0's result_y.pixels */
/*
	for(int i=0; i<num_ranks; ++i )
	{
		recv_counts[i] = (int *)((sobel_y.ysize/num_ranks) * sobel_y.xsize*sobel_y.depth) ;
		displs[i] = i*(sobel_y.ysize/num_ranks) * sobel_y.xsize*sobel_y.depth  ;
			
	}
*/	

	MPI_Gatherv(sobel_y.pixels,
		gather_size,
		MPI_CHAR,
		result_y.pixels,
		recv_counts,
		displs,
		MPI_CHAR,
		0,
		MPI_COMM_WORLD);

	convolve_time=MPI_Wtime();


	if (rank==0) {
//		store_jpeg("outx.jpg",&result_x);
//		store_jpeg("outy.jpg",&result_y);
	}


       /* combine to form output */
	if (rank == 0){
		
		combine(&result_x,&result_y,&new_image);
	        combine_time=MPI_Wtime();
		
        	/* Write data back out to disk */
		store_jpeg("out.jpg",&new_image);

	        store_time=MPI_Wtime();
	}

	/* Print stats */
	if (rank == 0){
		printf("Load time: %lf\n",load_time-start_time);
		printf("Convolve time: %lf\n",convolve_time-load_time);
		printf("Combine time: %lf\n",combine_time-convolve_time);
		printf("Store time: %lf\n",store_time-combine_time);
		printf("Total time = %lf\n",store_time-start_time);
	}

	MPI_Finalize();

	return 0;
}




