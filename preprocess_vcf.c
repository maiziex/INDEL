#include "common.h"

typedef struct
{
    int64_t candidate_num;
    int64_t max_candidate_num;
    int64_t* candiate_loci;
    int** GT;
} block_t;

typedef struct
{
    int64_t sample_num;
    int sample_name_max_length;
} pars_t;



block_t* init_block(int64_t max_candidate_num)
{
    block->candidate_num=0;
    block->max_candidate_num = max_candidate_num;
    block->candidate_loci = (int64_t*) malloc(block->max_candidate_num*sizeof(int64_t));
    block->GT = (int**) malloc(block->max_candidate_num*sizeof(int*));\
    for(i=0;i<block->max_candidate_num;i++)
    {
        block->GT[i]=(int*) malloc(pars->sample_num*sizeof(int));
    }
 
}

void recalloc_gl_block(block_t* block, pars_t* pars)
{
	block->candidate_loci = (int64_t*) realloc(block->candidate_loci, block->max_candidate_num * 2 * sizeof(int64_t));

	int64_t i, j;
	for(i = block->max_candidate_num; i < 2 * block->max_candidate_num; i++)
    {
            block->GT[i] = (int*) malloc(pars->sample_num * sizeof(int));
    }

	block->max_candidate_num *= 2;
	//printf("block size increase to %d\n", block->max_candidate_num);
}

void add_candidate_into_gl_block(block_t* block, pars_t* pars, int64_t locus,int** saved_GT)
{
	block->candidate_loci[block->candidate_num] = locus;
	int i;
	/*for(i = 0; i < pars->sample_num; i++)
	{
        block->GT[block->candidate_num][i] = saved_GT[block->candidate_num][i];
	}*/

	block->candidate_num++;
	if(block->candidate_num >= block->max_candidate_num)
	{
		recalloc_gl_block(block, pars);
	}
}



void init_query_samples(char* input_filename, pars_t* pars)
{
	FILE* input = fileOpenR(input_filename);

	int field_index = 0, sample_name_length;
        pars->sample_name_max_length = 0;
        char c, c_prev;

	while((c = (char) getc(input)) != EOF)
        {
                if(c == '#')
                {
                        c = (char) getc(input);
                        if(c == '#')
                        {
                                while(c != '\n') c = (char) getc(input);
                        }
                        else
                        {
                                while(c != '\n')
                                {
                                        c_prev = c;
                                        c = (char) getc(input);
                                        if(c == '\n') break;
                                        if((c_prev == ' ' || c_prev == '\t') && (c != ' ' && c != '\t'))
                                        {
                                                field_index++;
                                                if(field_index > VCF_MANDATORY)
                                                {
                                                        sample_name_length = 1;
                                                }
                                        }
                                        else if((c_prev != ' ' || c_prev != '\t') && (c != ' ' && c != '\t'))
                                        {
                                                if(field_index > VCF_MANDATORY)
                                                {
                                                        sample_name_length++;
                                                }
                                        }
                                        else
                                        {
                                                if(field_index > VCF_MANDATORY && sample_name_length > pars->sample_name_max_length)
                                                {
                                                        pars->sample_name_max_length = sample_name_length;
                                                }
                                        }
                                }
				break;
                        }
                }
	}

	pars->sample_num = field_index - VCF_MANDATORY;
	fileClose(input);
}





void read_vcf_data(CHAR *input_filename,pars_t* pars)
{
    init_query_samples(input_filename,pars);
    char c;
    block_t* block;
    init_block(100);
    FILE* input = fileOpenR(input_filename);
    while( (c= (char) getc(input)) != EOF)
    {    if(c == '#')
        {
            c = (char) getc(input);
            if(c == '#')
            {
                while(c!='\n') c = (char) getc(input);
            }
            else
            {
                while(c!='\n') c= (char) getc(input);
                break;
            }
        }
    }

    int i, j;
	int field_length; 
	int64_t locus;
    int max_candidate_num = 100;
    int num_of_semicolon_before_AF,numof_colon_before_GT,freeze; 
    char* chr = calloc(MAX_CHROMOSOME_LENGTH, sizeof(char));
	char* ref = calloc(MAX_CHROMOSOME_LENGTH, sizeof(char));
    char* alt = calloc(MAX_CHROMOSOME_LENGTH, sizeof(char));
    float* AF = calloc(MAX_LOCUS_LENGTH,sizeof(float));
    /*saved_GT = (int**) malloc(block->mm*sizeof(int*));
    for(i=0;i<block->max_candidate_num;i++)
    {
        saved_GT[i]=(int*) malloc(pars->sample_num*sizeof(int));
    }*/
 


    c = (char) getc(input);

	while(c != EOF)
	{
		memset(chr, 0, MAX_CHROMOSOME_LENGTH * sizeof(char));
		field_length = 0;
		while(c != '\t')
		{
			chr[field_length] = c;
			c = (char) getc(input);
			field_length++;
		}

		locus = 0;
		while((c = (char) getc(input)) != '\t')
		{
			locus = locus * 10 + c - '0';
		}

		//printf("%s:%d\n", chr, locus);
		while((c = (char) getc(input)) != '\t');
 	
        memset(ref, 0, MAX_CHROMOSOME_LENGTH * sizeof(char));
		field_length = 0;
		while((c = (char) getc(input)) != '\t')
		{
			ref[field_length] = c;
			field_length++;
		}

		memset(alt, 0, MAX_CHROMOSOME_LENGTH * sizeof(char));
		field_length = 0;
        while((c = (char) getc(input)) != '\t')
        {
            alt[field_length] = c;
            field_length++;
        }

		for(i = 0; i < 2; i++)
		{
			while((c = (char) getc(input)) != '\t');
		}

        
        freeze=0;
        num_of_semicolon_before_AF=0; 
		c = (char) getc(input);
		while(c != '\t')
		{
			c_prev = c;
			c = (char) getc(input);
			if(c_prev == ';' && !freeze) num_of_semicolon_before_AF++;
			if(c_prev == 'A' && c == 'F') freeze = 1;
            if(freeze==1)
            {
                c = (char) getc(input);
                AF[block->candidate_num]=(char) getc(input);
            }
		}

   	///////////////////////////////////////////////

        num_of_colon_before_GT = 0;
		c = (char) getc(input);
        while(c != '\t')
		{
			c_prev = c;
			c = (char) getc(input);
			if(c_prev == ':') num_of_colon_before_GT++;
			if(c_prev == 'G' && c == 'T') freeze = 1;
		}

		for(i = 0; i < pars->sample_num; i++)
		{
			for(j = 0; j < num_of_colon_before_GT; j++)
			{
				while(c != ':') c = (char) getc(input);
				c = (char) getc(input);
			}



            block->_GT[block->candidate_num][i]=0;
            while(c != ':' && c != '\t' && c != '\n')
            {
                c_prev = c;
                c = (char) getc(input);
                
                block->_GT[block->candidate_num][i] = block->GT[block->candidate_num][i] * 10 + c - '0';
                c = (char) getc(input);
            }
            c = (char) getc(input);

			///////////////////////////////////////////////

		}
        
        add_candidate_into_gl_block(block_t* block, pars_t* pars, int64_t locus,int** saved_GT)









     }
}


int main(int argc char *argv[])
{
    pars_t* pars =  (pars_t*) calloc(1,sizeof(pars_t));
   // set_default_pars(pars);
    read_vcf_data(argv[1], pars);

















}
:q

:q
