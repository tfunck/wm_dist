#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <pthread.h>
#include <stdbool.h>

#define TRUE 1
#define FALSE 0
#define MAX 99999 //FIXME: Only makes sense for brains at the mm scale, may have to be adapted for other applications
float pseudo_inf=999999; //9999999;
float dx, dy, dz, dx2, dy2, dz2;
int xmax, ymax, zmax, xymax;

static void pexit(char* string1, char* string2, int number) {
    fprintf(stderr, "Error: %s\t%s\n", string1, string2); 
    exit(number);
}

int real2voxel(float location, float min, float step){
//simple function to convert real to voxel coordinates
    //step=1;
    //if ( fmod((location-min), step) == 0 ) printf("Problem! %f\n",(location-min) / step );
    int out=(int) round((location-min)/step);
    //printf("(%f - %f) / %f --> %d\n", location,min, step, out);
    return(out );
}

typedef struct temp_point {
    int index;
	int id;
	int z;
	int y;
	int x;
	double wx;
	double wy;
	double wz;
	float value;

} point;


typedef struct image {
    char* filename;
    char* values_filename;
    //char* like_filename
    char** groups;
    double* data;
	point* points;
    unsigned int* count;
    double* start;
    double* step;
    double val;
    int n_roi;
    int ngroups;
    int ndim;
    int index;
    int n;
    int n3d;
    long int label;
    int tmax;
    int zmax;
    int ymax;
    int xmax;
    double xstep;
    double ystep;
    double zstep;
    double zstart;
    double ystart;
    double xstart;
    int tstep;
    int smooth;
} data;




struct node{
    float* dist;
    struct node** p;
    int index;
    //, z,y,x;
};


void swap(struct node* a, struct node* b){
    //int z=b->z;
    //int y=b->y;
    //int x=b->x;
    int index=b->index;
    
    float *dist=b->dist;
    struct node** p=b->p;
    //b
    //b->z=a->z;
    //b->y=a->y;
    //b->x=a->x;
    b->index=a->index;
    b->dist=a->dist;
    b->p=a->p;
    //a
    //a->z=z;
    //a->y=y;
    //a->x=x;
    a->index=index;
    a->dist=dist;
    a->p=p;
}

void sort_up(struct node* heap,int i, bool* index_array){
    int j=(int) floor(i/2);
    struct node* parent=&heap[j];
   
    if ( *parent->dist > *heap[i].dist && j !=0 ){   //parent is larger than child
    //printf("Sorting Up: Parent %d %f, Child %d %f\n", j,*parent->dist, j, *heap[i].dist  ); 
        index_array[parent->index]=i;
        index_array[heap[i].index]=j;
        swap(parent, &heap[i]);
        sort_up(parent,j, index_array);
    }
}
/*         15 (p)
//     	  /
//    	 27
//   	/  \
// (l) 25   30 (r)
*/


void sort_down(struct node* heap,int i, bool* index_array){
    int r=2*i;
    int l=2*i+1;
    int i2;
    struct node* rnode=&heap[r];
    struct node* lnode=&heap[l];
    struct node* temp;

    if ( *lnode->dist < *heap[i].dist || *rnode->dist < *heap[i].dist){
	    //if(i > max)	printf("Sorting Down: Current %d %f, LChild %d %f, RChild %d %f\n", i, *heap[i].dist, r, *rnode->dist, l, *lnode->dist  ); 
        if ( *lnode->dist < *rnode->dist ){   //lnode is smaller than right node
            temp=lnode;
            i2=l;
        }
        else {
            temp=rnode;
            i2=r;
        }
        //printf("\tReplacing %d,%f with %d,%f\n", i, *heap[i].dist, i2, *temp->dist);
        index_array[temp->index]=i;
        index_array[heap[i].index]=i2;
        swap(temp, &heap[i] );
        //printf("\tContinuig sort down at: %d %f\n", i2, *heap[i2].dist);
        sort_down(heap,i2, index_array);
    }

}

void sort(struct node* heap, int i, bool* index_array){
    int p=(int) floor(i/2);
    int l=i*2;
    int r=2*i+1;
    
    if(*heap[r].dist < *heap[i].dist || *heap[l].dist < *heap[r].dist ){
        sort_down(heap, i, index_array);
    }
    if( *heap[i].dist < *heap[p].dist){
        sort_up(heap, i, index_array);
    }

}

void delete(struct node* heap,int i, unsigned long n, bool* index_array, int run){
    //printf("Delete %d %d %f\n", heap[i].index, heap[i].dist, *heap[i].dist );
    index_array[heap[i].index]=n; //Set index of node in heap[i] to n (number of considered voxels
    index_array[heap[n].index]=i; //Set index of node heap[n] in index_array to i (1)
    swap(&heap[i],&heap[n]); //swap the nodes at position i and n
    heap[n].dist=&pseudo_inf; //set the distance of the node that is now at n to be equal to (pseudo) infinity. this means that this node has been travelled
	heap[n].p=NULL; //delete the parent node for the node at n
	//FIXME shouldn't the node be actually deleted from memory?
	//if(run >= 2285 && n >= 2858   ){ // 2298
	//  for(int o=0; o<n ; o++) 
	//	  if(heap[o].dist == NULL) printf("Uh oh, heap corrupted at %d\n", o);
    //}
    sort_down(heap, i, index_array); //sort the children of the first node

}

void insert(struct node* heap, int i, float* dist, int index, int z, int y, int x, bool* index_array){
    heap[i].index=index;
    //heap[i].z=z;
    //heap[i].y=y;
    //heap[i].x=x;
    heap[i].dist=dist;
    index_array[index]=i;
    //printf("Inserting %d %d %f\n", i, index, *heap[i].dist);
    sort_up(heap, i, index_array);
}
 
void* wm_dist_threaded(void*);
void useage();



struct wm_vol_args{
    char* example_fn;
    char* density_fn;
    data* img;
    unsigned int* density;
    int** gm_border;
    long int* img_vol;
    int write_vertex;
    long int label;
    int thread;
    int nthreads;
    unsigned long n;
    float* mat; 
};



int check_for_wm(long int* img_vol, int z, int y, int x, int zmax, int ymax, int xmax, int WM){
    for( int i=-1; i <= 1; i++ )
        for( int j=-1; j <= 1; j++ ) 
            for( int k=-1; k <= 1; k++ ){
                int zi=z+i;
                int yi=y+j;
                int xi=x+k; 
                if( img_vol[zi*ymax*xmax+yi*xmax+xi ] == WM) return(1);
            }
    return(0);
}

int** wm_gm_border(data* img, float** mesh, const long int label, const long int* img_vol,int** fill_wm, const unsigned long n, int* nReplace, int wm_search_depth){
    int ymax=img->ymax;
    int xmax=img->xmax;
    float zmin=img->zstart;
    float ymin=img->ystart;
    float xmin=img->xstart;
    float zstep=img->zstep;
    float ystep=img->ystep;
    float xstep=img->xstep;
    *nReplace=0;
    int** border=malloc(sizeof(*border)*n);
    /*****************
    *   Find border  * 
    ******************/
    for(int i=0; i<n; i++){
        border[i]=malloc(sizeof(**border)*5);
        float z0=mesh[i][0];
        float y0=mesh[i][1];
        float x0=mesh[i][2];
        int z=real2voxel(z0, zmin, zstep);
        int y=real2voxel(y0, ymin, ystep);
        int x=real2voxel(x0, xmin, xstep);
		//printf("z %f %f\n", zmin, zstep);
		//printf("y %f %f\n", ymin, zstep);
		//printf("x %f %f\n", xmin, zstep);
        int found_wm=FALSE;
        int p=0;
        float min_d=pseudo_inf;
        border[i][0]=z*xmax*ymax+y*xmax+x;
        border[i][1]=z;
        border[i][2]=y;
        border[i][3]=x;
        border[i][4]=1; //number of vertices with this start location
        while(found_wm==FALSE && p <= wm_search_depth){
            p++;
            for(int zi=z-p; zi < z+p; zi++ ){
                for(int yi=y-p; yi < y+p; yi++ ){
                    for(int xi=x-p; xi < x+p; xi++ ){
                        int index=zi*xymax+yi*xmax+xi;
                        int val=(int)img_vol[index];
						
                        if( (int) val== (int) label){
                            float d=(((float) zi)-z0)*(((float) zi)-z0) + (((float) yi)-y0)*(((float) yi)-y0)+ (((float) xi)-x0)*(((float) xi)-x0);
                            if(d < min_d){
                                min_d=d;
                                border[i][0]=index;
                                border[i][1]=zi;
                                border[i][2]=yi;
                                border[i][3]=xi;
                                border[i][4]=1; //Skip = FALSE
								if( zi==75 && yi == 81 && xi == 74) {
									printf("%d %d %d %d %d = %d\n", zi, xymax, yi, xmax, xi,zi*xymax+yi*xmax+xi );
									printf("%d %d %d %d --> %d\n", zi,yi,xi, index, img_vol[index]);
								}
                                found_wm=TRUE;   
                            }
                        } 
                    }
                }
            }
        }
    }
    for(int i=0; i<n; i++){ //For each voxels nearest to a vertex
        if(border[i][4] > 0){
            for(int j=0; j<n; j++){ //For all other voxels nearest to a vertex
                if(i != j && border[i][0]==border[j][0]){//If these two voxels are the same...
                    fill_wm[*nReplace][0]=i;//then keep track of the voxel we will copy from 
                    fill_wm[*nReplace][1]=j;//and the voxel we will copy to when filling in gaps in distance matrix.
                    border[i][4] += 1; //Increment the number of replicate voxels for this specific voxel
                    border[j][4]=0; //For the j voxel, set number of replicates to zero so that we know it's distance map is
                                    //already taken care of by voxel i
                    *nReplace += 1; //Increase the number of voxels to be replaced by 1/
                }
            }
        }
    }
    //for(int i=0; i<n; i++){ 
	//  	if( border[i][0] < 0 ) printf("%d\n", border[i][0]); 
	//}
	  //For each voxels nearest to a vertex
    //printf("From\tTo\n");
    //for(int i=0; i< *nReplace; i++){
    //    printf("%d\t%d\n", fill_wm[i][0], fill_wm[i][1]);
    //}
    //for(int i=0; i<n; i++) if(border[i][4] > 0) printf("%d\n", border[i][4]);
   printf("%3.1f% ( %d / %d ) of voxels can be skipped \n", 100.0 * (float) *nReplace / n, (int) *nReplace ,n); 

    return(border);
}

float** readVertices(const char* Meshfilename, const char* surface_mask_fn, unsigned long* nvertices, int subsample, int subsample_factor ){
    char *maskBuffer=NULL, *meshBuffer=NULL;
    int i=0; 
    int vertices;
    float** mesh;
    FILE *MeshFile=fopen(Meshfilename, "rt");
    FILE *MaskFile;
    char *token;
    char dlm[]=" ";
    int nmax, meshBufferN=0, maskBufferN=0; 
    int maskValue=1, recorded_vertices=0;

    getline(&meshBuffer, &meshBufferN, MeshFile) ;
    //read nvertices from file
    strtok(meshBuffer, dlm);  
    strtok(NULL, dlm); 
    strtok(NULL, dlm);  
    strtok(NULL, dlm); 
    strtok(NULL, dlm); 
    strtok(NULL, dlm); 
    *nvertices=atoi(strtok(NULL, dlm));
    
    if(subsample > 0) nmax=(int) 2 + round( (*nvertices-2)/pow(subsample_factor, subsample)); //2+(n1_vert-2)/(4^(s))
    else nmax=*nvertices ;

    if ( strcmp(surface_mask_fn, "") != 0  ) MaskFile=fopen(surface_mask_fn, "rt");

    if(nmax < 2 && subsample < 0){ 
        printf("Subsampled %d times with factor %d. This gives less than 2 vertices. Please try using less subsampling of the cortical surface mesh.\n"); 
        exit(1);
    }
    
    if(subsample > 0) printf("Mesh with %d vertices was subsampled %d (by factor of %d) to %d vertices\n", *nvertices, subsample, subsample_factor, nmax);
    mesh=malloc(sizeof(*mesh) * *nvertices);

    while( getline(&meshBuffer, &meshBufferN, MeshFile) != 0 )  {
        if ( strcmp(surface_mask_fn, "") != 0 ) { 
            getline(&maskBuffer, &maskBufferN, MaskFile);
            maskValue = atoi(maskBuffer);
        }

        if( maskValue != 0 ){
            mesh[recorded_vertices]=malloc(3*sizeof(**mesh));    
            mesh[recorded_vertices][2]=atof(strtok(meshBuffer, dlm)); //x
            mesh[recorded_vertices][1]=atof(strtok(NULL, dlm)); //y
            mesh[recorded_vertices][0]=atof(strtok(NULL, dlm)); //z
            recorded_vertices++;
        }

        i++;
        if( i >= nmax ) break; 
    }

    if(subsample > 0) *nvertices = nmax;

    if( strcmp(surface_mask_fn, "") != 0){ *nvertices=recorded_vertices;}

    free(maskBuffer); 
    free(meshBuffer); 
    fclose(MeshFile);
    return(mesh);
}


float quad(float a, float b, float c){
    float disc=b*b-4*a*c;
    float temp=(disc>0) ? sqrt(disc) : 0;
    float temp1=2*a;
    float sol1=(-b + temp ) / temp1;
    float sol2=(-b - temp ) / temp1;
    
	//printf("disc %f %f %f\n", disc, sol1, sol2);
	
	sol1 = (temp1 < 0 || disc < 0) ? 0 : sol1;
    sol2 = (temp1 < 0 || disc < 0) ? 0 : sol2;
    float sol=(sol1 > sol2) ? sol1 : sol2;
	return sol;
}

float quick_max(float a, float b, float c){
    if(b < a){
        if(a<c) return(c);
        else return(a);
    } else{
        if(b<c) return(c);
        else return(b);
    }
}

void update( bool* fixed, float* distances,  int index, int z, int  y, int  x){
  	int xpp, xmm, ypp, ymm, zpp, zmm;
    int ixp, ixpp, ixm, ixmm, iyp, iypp, iym, iymm, izp, izpp, izm, izmm;
    int xp, xm, yp, ym, zp, zm;
    float xdp, xdpp, xdm, xdmm,  ydp, ydpp, ydm, ydmm, zdm, zdmm, zdp, zdpp;
    float alist[6], blist[6], clist[6];
    float alta[3], altb[3], altc[3];
    float altdisc[3];
    float altsol[3];
    float sol, sol1, sol2;
    float alpha, disc, t, a,b, c; 
    float d; 
        d=distances[index]; 
		//printf("Distances: %d %d %d %f\n", z, y, x, d); //exit(0);
        xp=x+1;
        xpp=xp+1;
        xm=x-1;
        xmm=xm-1;
        yp=y+1;
        ypp=yp+1;
        ym=y-1;
        ymm=ym-1;
        zp=z+1;
        zpp=zp+1;
        zm=z-1;
        zmm=zm-1;
        ixp=z*xymax+y*xmax+xp; 
        ixpp=z*xymax+y*xmax+xpp;
        ixm=z*xymax+y*xmax+xm;
        ixmm=z*xymax+y*xmax+xmm;
        iyp=z*xymax+yp*xmax+x;
        iypp=z*xymax+ypp*xmax+x;
        iym=z*xymax+ym*xmax+x;
        iymm=z*xymax+ymm*xmax+x;
        izp=zp*xymax+y*xmax+x;
        izpp=zpp*xymax+y*xmax+x;
        izm=zm*xymax+y*xmax+x;
        izmm=zmm*xymax+y*xmax+x;
        xdp=distances[ixp];
        xdpp=distances[ixpp];
        xdm=distances[ixm];
        xdmm=distances[ixmm];
        ydp=distances[iyp];
        ydpp=distances[iypp];
        ydm=distances[iym];
        ydmm=distances[iymm];
        zdm=distances[izm];
        zdmm=distances[izmm];
        zdp=distances[izp];
        zdpp=distances[izpp];

        alist[0]=alist[1]=alist[2]=alist[3]=alist[4]=alist[5]=0;
        blist[0]=blist[1]=blist[2]=blist[3]=blist[4]=blist[5]=0;
        clist[0]=clist[1]=clist[2]=clist[3]=clist[4]=clist[5]=0;
        /*
        * Calculate new distance based on Eikonal equation
        */

        if( fixed[ixm]==TRUE) {//X-- second order
            if( fixed[ixmm]==TRUE){ 
                alpha = 9 / (4 * dx2);
                t=(4*xdm-xdmm)/3;
                alist[0] = alpha;
                blist[0] = -2*alpha*t;
                clist[0] = alpha*t*t;
				//printf("second order: %f %f\n", t, alpha);
            }
            else{ //X- first order
                alpha = 1/dx2;
                alist[0] = alpha;
                blist[0] = -2*xdm*alpha;
                clist[0] = xdm*xdm*alpha;
				//printf("first order\n");
            }
        } 
        if( fixed[ixp]==TRUE) {//X++ second order
            if( fixed[ixpp]==TRUE){ 
                alpha = 9 / (4 * dx2);
                t=(xdpp-4*xdp)/3;
                alist[1] = alpha;
                blist[1] = 2*alpha*t;
                clist[1] = alpha*t*t;
            }
            else{ //X+ first order
                alpha = 1/dx2;
                alist[1] = alpha;
                blist[1] = -2*xdp*alpha;
                clist[1] = xdp*xdp*alpha;
            }
        }
        if( fixed[iym]==TRUE) {//Y-- second order
            if( fixed[iymm]==TRUE){ 
                alpha = 9 / (4 * dy2);
                t=(4*ydm-ydmm)/3;
                alist[2] = alpha;
                blist[2] = -2*alpha*t;
                clist[2] = alpha*t*t;
            }
            else{ //Y- first order
                alpha = 1/dy2;
                alist[2] = alpha;
                blist[2] = -2*ydm*alpha;
                clist[2] = ydm*ydm*alpha;
            }
        } 
        if(fixed[iyp]==TRUE) {//Y++ second order
            if( fixed[iypp]==TRUE){ 
                alpha = 9 / (4 * dy2);
                t=(ydpp-4*ydp)/3;
                alist[3] = alpha;
                blist[3] = 2*alpha*t;
                clist[3] = alpha*t*t;
            }
            else{//Y+ first order
                alpha = 1/dy2;
                alist[3] = alpha;
                blist[3] = -2*ydp*alpha;
                clist[3] = ydp*ydp*alpha;
            }
        }        
        
        if(fixed[izm]==TRUE) {//Z-- second order
            if( fixed[izmm]==TRUE){ 
                alpha = 9 / (4 * dz2);
                t=(4*zdm-zdmm)/3;
                alist[4] = alpha;
                blist[4] = -2*alpha*t;
                clist[4] = alpha*t*t;
            }
            else{ //Z- first order
                alpha = 1/dz2;
                alist[4] = alpha;
                blist[4] = -2*zdm*alpha;
                clist[4] = zdm*zdm*alpha;
            }
        } 
        if(fixed[izp]==TRUE) {//Z++ second order
            if( fixed[izpp]==TRUE){
                alpha = 9 / (4 * dz2);
                t=(zdpp-4*zdp)/3;
                alist[5] = alpha;
                blist[5] = 2*alpha*t;
                clist[5] = alpha*t*t;
            }
            else{ //Z+ first order
                alpha = 1/dz2;
                alist[5] = alpha;
                blist[5] = -2*zdp*alpha;
                clist[5] = zdp*zdp*alpha;
            }
        }    
        /***************************
        *Try to find a 3D solution*
        ***************************/ 
        a=alist[0]+alist[1]+alist[2]+alist[3]+alist[4]+alist[5]; 
        b=blist[0]+blist[1]+blist[2]+blist[3]+blist[4]+blist[5];
        c=-1+clist[0]+clist[1]+clist[2]+clist[3]+clist[4]+clist[5];
		//printf("(a,b,c) = %f %f %f\n", a, b, c);
        sol=quad(a, b, c);
      	
		//printf("Sol 1 : %f\n", sol); 
        if( sol <= 0 || sol > pseudo_inf ){
            /***************************
             *Try to find a 2D solution*
             ***************************/
            alta[0]=alist[2]+alist[3]+alist[4]+alist[5];
            alta[1]=alist[0]+alist[1]+alist[4]+alist[5];
            alta[2]=alist[0]+alist[1]+alist[2]+alist[3];
            altb[0]=blist[2]+blist[3]+blist[4]+blist[5];
            altb[1]=blist[0]+blist[1]+blist[4]+blist[5];
            altb[2]=blist[0]+blist[1]+blist[2]+blist[3];
            altc[0]=-1+clist[2]+clist[3]+clist[4]+clist[5];
            altc[1]=-1+clist[0]+clist[1]+clist[4]+clist[5];
            altc[2]=-1+clist[0]+clist[1]+clist[2]+clist[3];

            altsol[0] = quad(alta[0], altb[0], altc[0]);
            altsol[1] = quad(alta[1], altb[1], altc[1]);
            altsol[2] = quad(alta[2], altb[2], altc[2]);

            /*if(x==111 && y==159 && z==82){ 
                printf("2D Solutions\n");
                for(int kk=0; kk<3; kk++) printf("\t%f %f %f\n", alta[kk], altb[kk], altc[kk]);
                printf("%f  %f  %f \n", altsol[0],  altsol[1], altsol[2]); 
            }*/

            sol=quick_max(altsol[0], altsol[1], altsol[2] ); 
			//printf("Sol 2 : %f\n", sol); 
            if (  sol <=0 || sol > pseudo_inf  ){
                /***************************
                 *Try to find a 1D solution*
                 ***************************/ 
                alta[0]=alist[2]+alist[3];
                alta[1]=alist[0]+alist[1];
                alta[2]=alist[4]+alist[5];
                altb[0]=blist[2]+blist[3];
                altb[1]=blist[0]+blist[1];
                altb[2]=blist[4]+blist[5];
                altc[0]=-1+clist[2]+clist[3];
                altc[1]=-1+clist[0]+clist[1];
                altc[2]=-1+clist[4]+clist[5];
                
                altsol[0] = quad(alta[0], altb[0], altc[0]);
                altsol[1] = quad(alta[1], altb[1], altc[1]);
                altsol[2] = quad(alta[2], altb[2], altc[2]);
                sol=quick_max(altsol[0], altsol[1], altsol[2]);
				
				//printf("Sol 3: %f\n", sol );	
            }
        } 
		//printf("Final Sol: %f\n", sol);
		if( isnan(sol) || sol > pseudo_inf || sol <= 0  ) { 
            //printf("Broken sol: %f %f\n", sol,distances[index] ); 
            //exit(0); 
        }
        else {
		 	distances[index]= sol ;
		}
        //distances[index]= (sol > 0) ? sol : distances[index];
}
int add_neighbours(struct node* considered, long int* img_vol, float* distances, bool* fixed, bool* considered_array, int* nconsidered, const long int label,const int  z, const int  y, const int x, const int init ){
    int i0, i1, i2, i3, i4, i5;
    int    xp=x+1;
    int    xm=x-1;
    int    yp=y+1;
    int    ym=y-1;
    int    zp=z+1;
    int    zm=z-1;
	i0=z*xymax+y*xmax+xp;
	i1=z*xymax+y*xmax+xm;
	i2=z*xymax+yp*xmax+x;
	i3=z*xymax+ym*xmax+x;
	i4=zp*xymax+y*xmax+x;
	i5=zm*xymax+y*xmax+x;
	
	/*printf("\nInitialization: %d\n", init);

	  printf("%d %d %d\n", zmax, ymax, xmax);
	  printf("%d %d %d %d %d %d %f\n",z,y, x, img_vol[z*xymax+y*xmax+x], fixed[z*xymax+y*xmax+x], considered_array[z*xymax+y*xmax+x], distances[z*xymax+y*xmax+x] );
	  printf("%d %d %d %d %d %d %f\n",z,y,xp, img_vol[i0], fixed[i0], considered_array[i0], distances[i0]);
	  printf("%d %d %d %d %d %d %f\n",z,y,xm, img_vol[i1], fixed[i1], considered_array[i1], distances[i1]);
	  printf("%d %d %d %d %d %d %f\n",z,yp,x, img_vol[i2], fixed[i2], considered_array[i2], distances[i2]);
	  printf("%d %d %d %d %d %d %f\n",z,ym,x, img_vol[i3], fixed[i3], considered_array[i3], distances[i3]);
	  printf("%d %d %d %d %d %d %f\n",zp,y,x, img_vol[i4], fixed[i4], considered_array[i4], distances[i4]);
	  printf("%d %d %d %d %d %d %f\n",zm,y,x, img_vol[i5], fixed[i5], considered_array[i5], distances[i5]);
	  */
	if ( (int) img_vol[i0]==(int)label && fixed[i0]==FALSE  ){
		if(considered_array[i0] == FALSE) {
			*nconsidered += 1;
			 insert(considered, *nconsidered, &(distances[i0]), i0, z, y, xp, considered_array );
		}
		if(init) distances[i0]=dx;
		else{  
			update(fixed,distances, i0, z, y, xp);
			sort(considered, considered_array[i0], considered_array);
		}
	}
	if ((int)img_vol[i1]==(int)label && fixed[i1]==FALSE  ){ 
		if(considered_array[i1] == FALSE) {
			*nconsidered += 1;
			insert(considered,  *nconsidered, &(distances[i1]), i1, z, y, xm, considered_array );
		}
		if(init) distances[i1]=dx;
		else{ 
			update(fixed,distances, i1,z , y,  xm);
			sort(considered, considered_array[i1], considered_array);
		}
	}

	if ((int)img_vol[i2]==(int)label && fixed[i2]==FALSE   ){
		if(considered_array[i2] == FALSE) {    
			*nconsidered += 1;
			insert(considered,  *nconsidered, &(distances[i2]), i2, z, yp, x, considered_array );
		}
		if(init) distances[i2]=dy;
		else{ 
			update(fixed,distances, i2, z, yp, x);
			sort(considered, considered_array[i2], considered_array);
		}
	}
	if ((int)img_vol[i3]==(int)label && fixed[i3]==FALSE   ){ 
		if(considered_array[i3] == FALSE) {    
			*nconsidered += 1;
			insert(considered,  *nconsidered, &(distances[i3]), i3, z, ym, x, considered_array);
		}
		if(init) distances[i3]=dy;
		else{ 
			update(fixed,distances, i3,z , ym, x);
			sort(considered, considered_array[i3], considered_array);
		}
	}

	if ((int)img_vol[i4]==(int)label && fixed[i4]==FALSE  ){ 
		if(considered_array[i4] == FALSE) {
			*nconsidered += 1;
			insert(considered,  *nconsidered, &(distances[i4]), i4, zp, y, x, considered_array);     

		}
		if(init) distances[i4]=dz;
		else{ 
			update(fixed,distances, i4, zp , y, x);
			sort(considered, considered_array[i4], considered_array);
		}
	}

	if ((int)img_vol[i5]==(int)label && fixed[i5]==FALSE  ){
		if(considered_array[i5] == FALSE) {    
			*nconsidered += 1;
			insert(considered,  *nconsidered, &(distances[i5]), i5, zm, y, x, considered_array); 
		}
		if(init) distances[i5]=dz;
		else{ 
			update(fixed,distances, i5, zm , y, x);
			sort(considered, considered_array[i5], considered_array);
		}
	}
	/*
	  printf("%d %d %d %d %d %d %f\n",z,y, x, img_vol[z*xymax+y*xmax+x], fixed[z*xymax+y*xmax+x], considered_array[z*xymax+y*xmax+x], distances[z*xymax+y*xmax+x] );
	  printf("%d %d %d %d %d %d %f\n",z,y,xp, img_vol[i0], fixed[i0], considered_array[i0], distances[i0]);
	  printf("%d %d %d %d %d %d %f\n",z,y,xm, img_vol[i1], fixed[i1], considered_array[i1], distances[i1]);
	  printf("%d %d %d %d %d %d %f\n",z,yp,x, img_vol[i2], fixed[i2], considered_array[i2], distances[i2]);
	  printf("%d %d %d %d %d %d %f\n",z,ym,x, img_vol[i3], fixed[i3], considered_array[i3], distances[i3]);
	  printf("%d %d %d %d %d %d %f\n",zp,y,x, img_vol[i4], fixed[i4], considered_array[i4], distances[i4]);
	  printf("%d %d %d %d %d %d %f\n",zm,y,x, img_vol[i5], fixed[i5], considered_array[i5], distances[i5]);
	  */
    return(*nconsidered);
}

int min_path(int index, int z, int y, int x, 
        float* distances /*distance map*/, 
        unsigned int* density, /*map of voxels tranversed by min distance paths*/ 
        _Bool* fixed, /*map of fixed points*/ 
        long int* img_vol, 
        float min, 
        int factor){
    int min_index=pseudo_inf;
    float local_min=min;
    int min_zi, min_yi, min_xi;
    int break_point=0, c=1;
    

    //while(break_point == 0 || c < 3 ){
        int z0 = (z-c >= 0) ? z-1 : 0;
        int y0 = (y-c >= 0) ? y-1 : 0;
        int x0 = (x-c >= 0) ? x-1 : 0;
        int z1 = (z+c < zmax) ? z+1 : zmax;
        int y1 = (y+c < ymax) ? y+1 : ymax;
        int x1 = (x+c < xmax) ? x+1 : xmax;
        for(int zi=z0; zi <= z1; zi++){ 
            for(int yi=y0; yi <= y1; yi++){ 
                for(int xi=x0; xi <= x1; xi++){
                    int index1=zi*xymax+yi*xmax+xi;
                    if( fixed[index1] == TRUE && index1 != index && distances[index1] < local_min){
                        min_index=index1;
                        local_min=distances[index1]; 
                        min_zi=zi;
                        min_yi=yi;
                        min_xi=xi;
                        break_point=1;
                    }
                } 
            }
        }
     //   c++;
    //}
    
    if(min_index != pseudo_inf && break_point == 1 ) {
        //printf("\t%d %d %d : %d %d %f %d, %d\n", min_zi, min_yi, min_xi, fixed[min_index], distances[min_index], img_vol[min_index], density[min_index], c);
        density[min_index] += factor;
        min_path(min_index, min_zi, min_yi, min_xi, distances, density, fixed, img_vol, local_min, factor);
   }

    return(0);
}

int min_paths(int cur_i, float* distances, unsigned int* density, _Bool* fixed, long int* img_vol, int** gm_border, unsigned long n, int factor){
    //FIXME: Seems like the algorithm hits a local minimum? 
    
    for(int i=0; i<n; i++){
        //For each voxel that is on the GM-WM border, track the minimum distance back to the curent
       //voxel on the GM-WM border. 
        if(gm_border[i][4] != TRUE && i != cur_i){
            int index=gm_border[i][0];
            int z=gm_border[i][1];
            int y=gm_border[i][2];
            int x=gm_border[i][3];
            float min=distances[index];
            min_path(index, z, y, x, distances, density, fixed, img_vol, min , factor);
        } /**/ //else printf("Skipped\n");
    }
    return(0);
}

int eikonal(struct node* considered, int*  img_vol, float*  distances, bool* fixed, bool* considered_array, const long int label, int  z, int y, int  x, int run){
    int nconsidered=0; 
    /************************************
     * Initialize the considered points *
     ************************************/
	
    add_neighbours(considered, img_vol, distances,fixed, considered_array, &nconsidered, label,  z, y, x, 1);
	//printf("Distances A: %f\n", distances[z*xymax+y*xmax+x]);

	while(nconsidered > 0){
        /********************************************************
         * 3. Add mininum distance point to list of fixed points
         ********************************************************/
        int index=considered[1].index;
        z=(int) floor(index / (xymax));
        y=(int) floor(index-  z*xymax)/xmax; 
        x=(int) floor(index-z*xymax-y*xmax);
       	//printf("%d %d %d -> %f\n", z, y, x, distances[index]);
		delete(considered, 1, nconsidered, considered_array, run);
        //Decrease the number of considered points
        nconsidered--;
        fixed[index]=TRUE;
        /***********************************************************************
         * Add points surrounding new fixed point to list of considered points
         ***********************************************************************/
        add_neighbours(considered, img_vol, distances,fixed, considered_array, &nconsidered,label,  z, y, x, 0);
    
		//printf("\n%d\n", nconsidered);
	}


    return(0);
}

void wm_dist_singlethread(data* img, long int* img_vol, int** gm_border, float* mat, const unsigned long n,const long int label, int write_vertex,char* example_fn, unsigned int *density, char* density_fn,  const  int start,const int step){
    int max =img->n3d;
    bool* fixed=malloc(max * sizeof(*fixed));
    bool* considered_array=malloc(max * sizeof(*considered_array));
    struct node* considered=malloc(max *sizeof(*considered)); //, *considered_last=FALSE, *minNode=FALSE;
    float* distances=malloc(max*sizeof(*distances));
    int index;
    int x, y, z;
    int i;

    for( i=start; i < n; i += step){ //Iterate over nodes on WM-GM border
		printf("\rThread %d: %3.1f",start, (float) 100.0 * i/n); fflush(stdout);
        for(int j=0; j<max; j++){ 
            considered[j].dist=&pseudo_inf;
            distances[j]=pseudo_inf;
            fixed[j]=considered_array[j]=FALSE;
        }
        if(gm_border[i][4]==0){ continue;}
        index=gm_border[i][0];
        fixed[index]=TRUE; //Set the starting fixed node
        z=gm_border[i][1];
        y=gm_border[i][2];
        x=gm_border[i][3];
        distances[index]=0.0;

		eikonal(considered, img_vol, distances,fixed, considered_array, label,  z, y, x, i);
		
		if( strcmp(density_fn,  "") != 0 ){
            // Find minimum path across distance gradient 
			min_paths(i, distances, density, fixed, img_vol, gm_border, n, gm_border[i][4]);
        }

        if(  (i ==write_vertex || write_vertex == -1) && strcmp(example_fn, "") != 0 ){
            for(int j=0; j<max; j++) distances[j] *= fixed[j];
            if(write_vertex==-1) sprintf(example_fn, "vertex_%d.mnc", i);
            //printf("Writing vertex %d to %s\n", i, example_fn);
            //writeVolume(example_fn, distances, img->start, img->step, img->wcount, MI_TYPE_FLOAT );
        }
		if( mat != NULL  ){
            for(int v=0; v<n; v++){//for v in list of vertices
                z=gm_border[v][1];//get the index value of the initial wm voxel
                y=gm_border[v][2];
                x=gm_border[v][3];
                index=gm_border[v][0];
                mat[i*n+v]=distances[index];
            }
        }
    }
    free(considered);
    free(distances);
    free(fixed); 
    free(considered_array); 
}

void* wm_dist_threaded(void* args){
    data* img= ((struct wm_vol_args*) args)->img;
    long int label= ((struct wm_vol_args*) args)->label;
    long int* img_vol=((struct wm_vol_args*) args)->img_vol;
    float* mat=((struct wm_vol_args*) args)->mat;
    int** gm_border=((struct wm_vol_args*) args)->gm_border;
    unsigned long n=((struct wm_vol_args*) args)->n;
    int write_vertex=((struct wm_vol_args*) args)->write_vertex;
    char* example_fn=((struct wm_vol_args*) args)->example_fn;
    char* density_fn=((struct wm_vol_args*) args)->density_fn;
    unsigned int* density=((struct wm_vol_args*) args)->density;
    int thread= ((struct wm_vol_args*) args)->thread;
    int nthreads= ((struct wm_vol_args*) args)->nthreads;
    wm_dist_singlethread(img, img_vol, gm_border, mat, n, label, write_vertex, example_fn, density, density_fn, thread, nthreads);
}

int wm_dist_multithreaded(data* img, long int* img_vol, int** gm_border, long int label, float* mat, unsigned long n, int nthreads, int  write_vertex, char* example_fn, unsigned int* density, char* density_fn){
    int rc;

    pthread_t threads[nthreads];
    struct wm_vol_args thread_args[nthreads];
    for(int t=0; t< nthreads; t++){
        thread_args[t].img_vol=img_vol;
        thread_args[t].img=img;
        thread_args[t].label=label;
        thread_args[t].mat=mat;
        thread_args[t].n=n;
        thread_args[t].gm_border=gm_border;
        thread_args[t].write_vertex=write_vertex;
        thread_args[t].example_fn=example_fn;
        thread_args[t].density_fn=density_fn;
        thread_args[t].density=density;
        thread_args[t].thread=t;
        thread_args[t].nthreads=nthreads;
        rc= pthread_create(&threads[t], NULL, wm_dist_threaded, (void *) &thread_args[t] ); 
        if(rc!=0) pexit("Error: creating thread","" , 1);
    }
    for(int t=0; t<nthreads;t++){
        rc = pthread_join(threads[t], NULL);
        if(rc!=0) pexit("Error: joining thread", "", 1);
    }

    return(0);
}



/***************************************************************************************************
 *Author: Thomas Funck 
 *Contact: thomas.funck@mail.mcgill.ca
 *Date: December 1, 2016
 *Location: Montreal Neurological Institute, Montreal QC
 *
 *Name: wm_dist  
 *Synopsis: - program to calculate minimum distance from a region to any point within anatomically constrained space
 *
 *Description:
 *  Calculates distances from 3D region within anatomically constrained region. 
 *
 *To do:
 *  -  Nodes are currently only defined as closest voxels to vertices on mesh. 
 *      Should extend to add option for volumetric labels around which to create nodes.
 *      This would add possibility of measuring distances to/from subcortical gm
 *  - Calculate distances along cortical surface mesh for geodesic distances
 *  - Read/write nifti images for greater compatibility
 *  - Write python script
 *  - Python wrapper
 *
 *
 * ****************************************************************/
int wm_dist(long int* img_vol, char* mesh_fn,  char* matrix_fn, char* surface_mask_fn, char* example_fn, char* density_fn, long int label, unsigned int wm_search_depth, unsigned int max_threads, unsigned int write_vertex, unsigned int subsample, unsigned int subsample_factor , double* steps, long int* maxima, double* starts, int VERBOSE  ){
  	if (VERBOSE >= 2) {
	  printf("%d\n", img_vol); 
	  printf("%s\n", mesh_fn);   
	  printf("%s\n", matrix_fn); 
	  printf("%s\n", surface_mask_fn); 
	  printf("%s\n", example_fn); 
	  printf("%s\n", density_fn); 
	  printf("%d\n", label); 
	  printf("%d\n", wm_search_depth); 
	  printf("%d\n", max_threads); 
	  printf("%d\n", write_vertex); 
	  printf("Subsample %d\n", subsample); 
	  printf("Subsample factor %d\n", subsample_factor); 
	  printf("%f %f %f\n", steps[0], steps[1], steps[2]); 
	  printf("%d %d %d\n", maxima[0], maxima[1], maxima[2]);
	  printf("%f %f %f\n", starts[0], starts[1], starts[2]);  
	  printf("%d\n", VERBOSE);
	  printf("Matrix : %s\n", (char*) matrix_fn);	
	}
    float* mat=NULL;
    FILE* matrix_file=fopen(matrix_fn, "w+");
    data img;
    double* dist_vol;
    int** fill_wm;
    char *file_inputs[]={mesh_fn}; //{img.filename,  mesh_fn}; //={mesh_filename, node_values_filename};
    int n_input_files=2;
    int nthreads=sysconf(_SC_NPROCESSORS_ONLN);
    unsigned long n, nReplace=0;

    printf("Mesh:  \t%s\n", mesh_fn);
    printf("Label: \t%d\n", label);
    printf("Search depth: %d\n",wm_search_depth);
    if( strcmp(surface_mask_fn, "") != 0 )  printf("Surface Mask: %s\n", surface_mask_fn);  
    if( strcmp(example_fn, "") != 0 ) printf("Example_fn: %s\n", example_fn);  
    if( strcmp(density_fn, "") != 0 ) printf("Density_fn: %s\n", density_fn);  
    
    //Set maximum number of threads to use. 
    //Default is the number of cores on system
    if(max_threads < nthreads && max_threads > 0) nthreads=max_threads;
    if(VERBOSE) printf("Number of threads: %d\n", nthreads);


    float** mesh=readVertices(mesh_fn, surface_mask_fn, &n, subsample, subsample_factor);

    if(VERBOSE) printf("Number of nodes: %d\n", n);
	
	img.n3d = maxima[0]*maxima[1]*maxima[2];
    unsigned int* density=calloc(img.n3d, sizeof(*density));
    img.zstep = steps[0];
	img.ystep = steps[1];
	img.xstep = steps[2];
	
	dx=fabsf(img.xstep);
    dx2=dx*dx;
    dy=fabsf(img.ystep);
    dy2=dy*dy;
    dz=fabsf(img.zstep);
    dz2=dz*dz;

	img.zstart = starts[0];
	img.ystart = starts[1];
	img.xstart = starts[2];
	zmax = img.zmax = maxima[0];
    ymax = img.ymax = maxima[1];
    xmax = img.xmax = maxima[2];
    xymax=img.xmax*img.ymax;

    if(VERBOSE) printf("dx,dy,dz: %f %f %f\n", dx,dy,dz); 
    if( strcmp(matrix_fn, "") != 0) mat=malloc(sizeof(*mat) * n * n);
    else mat=NULL;
    
	fill_wm=malloc(sizeof(*fill_wm) *n);
    for(int i=0; i<n; i++) fill_wm[i]=malloc(sizeof(**fill_wm) *2);
	int * img_vol_0 = malloc( img.n3d *  sizeof(int)  );
	
	for(int z=0; z<img.zmax; z++ ){
	  	for(int y=0; y<img.ymax; y++ ){
            for(int x=0; x<img.xmax; x++ ){
                if(z<=1 || y <= 1 || x<=1 || z >= img.zmax-2 || y >= img.ymax-2 || x >= img.xmax-2){
				  	img_vol[z*img.ymax*img.xmax+y*img.xmax+x]=0;
                }
			}
		}
	}
    
	//if(matrix_fn != NULL) free(mat);

    /*Find GM-WM Border*/
    int** gm_border=wm_gm_border(&img, mesh,label, img_vol, fill_wm,  n, &nReplace, wm_search_depth);

		
	if(nthreads > 1){
        wm_dist_multithreaded(&img, img_vol, gm_border, label, mat, n, nthreads, write_vertex, example_fn, density, density_fn);
    } else {
        wm_dist_singlethread(&img, img_vol, gm_border, mat, n, label, write_vertex, example_fn, density, density_fn, 0, 1);
    }
    //if(density_fn != NULL){ 
        //int nkept=n-nReplace;
        //for(int i=0; i< img.n3d; i++) density[i] = ((float) density[i]/nkept);
    //    if(VERBOSE) printf("Writing density map to %s.\n", density_fn);
    //    writeVolume(density_fn, density, img.start, img.step, img.wcount, MI_TYPE_INT );
    //}

    if( strcmp(matrix_fn, "") != 0 ){
        for(int i=0; i<nReplace; i++){
            int from=fill_wm[i][0];
            int to=fill_wm[i][2];
            for(int j=0; j<n; j++){
                //printf("%d %f --> %f\n",i, mat[from*n+j],  mat[to*n+j] );
                mat[to*n+j]=mat[from*n+j];
            }
        }
    }

    /******************************
     * Write matrix to outputfile
     ******************************/
    if( strcmp(matrix_fn, "") != 0 ){
        if(VERBOSE) printf("\nDistances calculated.\nWriting to %s\n", matrix_fn);
        for(int i=0; i<n; i++){
            for(int j=0; j<n; j++){
			  	fprintf(matrix_file, "%f", mat[i*n+j]);
                if(j<n-1) fprintf(matrix_file, ",");
            }
            fprintf(matrix_file, "\n");
        }
    }
    for(int i=0; i<n; i++){ 
        free(fill_wm[i]);
        free(gm_border[i]);
        free(mesh[i]);
    }
    free(mesh);
    free(gm_border);
    free(fill_wm); 
    free(img_vol);
    free(density);
    return 0;
}




void useage(){
    printf("\nName:\n\twm_dist ~ calculate pairwise distances between points in gm through the  wm.\n");
    printf("Description:\n\tSolves the Eikonal equation for wavefront propagation using the Fast Marching algorithm\n");
    printf("Options:\n");
    printf("   Outputs:\n");
    printf("\t-density <.mnc>\t\t\t\tMINC volume with number of min paths through voxel, normalized to number of vertices\n");
    printf("\t-matrix <.csv/.txt>\t\t\tCSV file with pairwise vertex distances\n");
    printf("\t-write_vertex <int> <example.mnc>\tWrite distance for a vertex to MINC file.\n");
    printf("   Parameters:\n");
    printf("\t-wm_search_depth <int>\t\t\tSearch depth of voxels around a node to be checked for WM. Default=2\n");
    printf("\t-threads <int>\t\t\t\tNumber of threads for parrallel processing. Default=Number of cores\n");
    printf("\t-subsample <int>\t\t\tNumber of times to subsample cortical surface mesh.\n");
    printf("\t-subsample_factor <int>\t\t\tFactor by which to subsample (default=4 for CIVET obj).\n");
    printf("   Miscallaneous:\n");
    printf("\t-silent\t\t\t\t\tRun silently.\n");
    printf("Useage:\n\t<options> mask.mnc vertices.obj  <WM label> \n");
    printf("\tNote: options must come before mandatory arguments\n");
    printf("Example:\n\t1. -matrix wm_dist.csv sub01_pve_classify.mnc  wm_mesh.obj 3\n");
    printf("\t2. wm_dist -density density.mnc -threads 4 -write_vertex 105 distances.mnc wm_mask.mnc inner.obj 1\n");
    exit(1);
}
