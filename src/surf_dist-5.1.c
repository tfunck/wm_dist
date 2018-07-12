#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "minc2.h"
#define TRUE 1
#define FALSE 0
#define MAX 99999 //FIXME: Only makes sense for brains at the mm scale, may have to be adapted for other applications
int VERBOSE=FALSE;



struct temp_vertex;
struct temp_edge;

typedef struct temp_vertex {
  double z;
  double y;
  double x;
  double dist;
  double val;
  double previous_distance;
  //for voxels the pointer below should 
  //point to the nearest vertex. For vertices
  //the pointer should go to the ngh 
  //of that vertex
  struct temp_vertex *previous;
  //the info below only applies to vertices
  //on a mesh, not to voxels
  struct temp_vertex** ngh;
  struct temp_edge** incoming_edges;
  struct temp_edge** outgoing_edges;
  int n_incoming_edges;
  int n_outgoing_edges;
  int inside;
  int border;
  int nngh;
  int index;
  int n;
} vertex;


typedef struct temp_edge {
    vertex* start;
    vertex* end;
    double length;
    int travelled;

} edge;

double dist(double a, double b);
int check_if_border(vertex* c_vtx);
void useage();

void pexit(char* string,  char* file, int line, int number) {
    printf("%s in file %s at line %d.\n", string, line, number); 
    exit(number);
}

int length(char *array[]);
int interpolate_sphere( char * input_object_filename, int** n_ngh, int*** ngh ); 

vertex** read_vertices(char* mesh_filename, char* vertex_values_filename, int* nvertices, float inside_label){
    char buffer1[100],  buffer2[100]; //should probably improve buffer length so that it is not a fixed size
    int i; 
    FILE *mesh_file=fopen(mesh_filename, "rt");
    FILE *vertex_values_file=fopen(vertex_values_filename, "rt");
    char dlm[]=" "; 
    vertex* mesh;
    vertex** mesh_ptr;


    if(vertex_values_file==NULL) {
        printf("Error opening %s\n", mesh_filename);
        exit(1);
    }
    if(mesh_file==NULL) {
        printf("Error opening %s\n", mesh_filename);
        exit(1);
    }

    fgets(buffer1, sizeof(buffer1), mesh_file) ;
    //read nvertices from file
    strtok(buffer1, dlm);  
    strtok(NULL, dlm); 
    strtok(NULL, dlm);  
    strtok(NULL, dlm); 
    strtok(NULL, dlm); 
    strtok(NULL, dlm); 

    *nvertices=atoi(strtok(NULL, dlm));
    mesh=malloc(sizeof(vertex) * *nvertices);
    mesh_ptr=malloc(sizeof(vertex*) * *nvertices);
    for(i=0; i< *nvertices; i++){
        fgets(buffer1, sizeof(buffer1), mesh_file);
        fgets(buffer2, sizeof(buffer2), vertex_values_file); 
            mesh[i].x=atof(strtok(buffer1, dlm)); 
            mesh[i].y=atof(strtok(NULL, dlm)); 
            mesh[i].z=atof(strtok(NULL, dlm)); 
            mesh[i].index=i;
            mesh[i].val=atof(buffer2);
            mesh[i].dist=0;
            mesh[i].n_incoming_edges=mesh[i].n_incoming_edges=0;
            mesh[i].incoming_edges=mesh[i].outgoing_edges=NULL;
            mesh[i].border=FALSE;
            mesh[i].n=0;
            if(mesh[i].val == inside_label){
                mesh[i].inside=TRUE;
                mesh[i].dist=0;
            }
            else {
                mesh[i].dist=0;
                mesh[i].inside=FALSE;
            }
            //printf("%s\t%f\t", buffer2, mesh[i].val);
            mesh_ptr[i]=&mesh[i];
        }

    fclose(mesh_file);
    fclose(vertex_values_file);
    return(mesh_ptr);
}


int linkVerticesToNeighbours(vertex** mesh, int* n_ngh, int** ngh, int nvertices, int max_vertex_length){
    int i, k;
    int ngh_count;
    edge* e;
    vertex* c_vtx, *t_vtx;
    int j=0;
    for (i=0; i< nvertices; i++){
        c_vtx=mesh[i];
 	    c_vtx->ngh=malloc(sizeof(vertex*) * (n_ngh[i]) );
	    c_vtx->nngh=n_ngh[i];
	    c_vtx->ngh[0]=mesh[i];
	    //the first neighbour points to itself. this may seem strange but its easier to work with an array that contains all
	    //the points we are interested in.

	    for(k=0; k< n_ngh[i]; k++){
		    //because k starts at one we have to subtract k by 1 in ngh, otherwise we spill over into
		    //adjascent arrays
            //first check that the neighbour actually exists
            if(ngh[i][k] < nvertices  ){
		        c_vtx->ngh[k]=mesh[ ngh[i][k] ];
                t_vtx=c_vtx->ngh[k];
                
                //if we find a neighbouring vertex that is outside
                if(t_vtx->inside != TRUE){
                    double length=sqrt( dist(c_vtx->x, t_vtx->x) + dist(c_vtx->y, t_vtx->y) + dist(c_vtx->z, t_vtx->z)  );
                    if (length < max_vertex_length) {
                    //if the vertex mesh[i] is connected to a neigbour that is outside
                    //then create an edge connecting the two
                        e=malloc(sizeof(edge));
                        e->start=c_vtx;
                        e->end=t_vtx;
                        e->travelled=FALSE;
                        e->length=length;
                        c_vtx->outgoing_edges=realloc(c_vtx->outgoing_edges, sizeof(edge*) * (1+c_vtx->n_outgoing_edges));
                        c_vtx->outgoing_edges[c_vtx->n_outgoing_edges]=e;
                        c_vtx->n_outgoing_edges += 1;
                        t_vtx->incoming_edges=realloc(t_vtx->incoming_edges, sizeof(edge*) * (1+t_vtx->n_incoming_edges));
                        t_vtx->incoming_edges[t_vtx->n_incoming_edges]=e;
                        t_vtx->n_incoming_edges += 1;

                    } 
                    //set current vertex's border flag to true
                    if(c_vtx->inside == TRUE){ 
                        c_vtx->border=TRUE;
                        //printf("In: %d\tBorder:\t%d\n", c_vtx->inside, c_vtx->border);
                    }
                }
                //increment number of outgoing and incoming edges
       
            }
	    }
    }
    for (i=0; i< nvertices; i++) mesh[i]->border=check_if_border(mesh[i]);


    return(0);
}

int length(char *array[]) {
    int i;
    for(i=0; array[i] != NULL; i++);
    return(i);
}


int check_input_files(char *file_inputs[]){
    int n_file_inputs=length(file_inputs);
    int i;
    if(VERBOSE) printf("Number of file inputs: %d\n", n_file_inputs);
    for(i=0; i < n_file_inputs; i++){
        if( access( file_inputs[i], R_OK ) != -1 ) {
            if(VERBOSE) printf("Can read file: %s\n", file_inputs[i]);
        } else {
            printf("Error: could not access %s\n", file_inputs[i]);
            return 1;
        }
    }
    return 0 ;
}

double dist(double a, double b){
    return ((a-b)*(a-b));
}


int find_euclidian_distance(vertex** mesh, int nvertices, int inside_label){
    double x_sum=0, y_sum=0, z_sum=0;
    double euclidian_distance;
    int n=0;
    vertex* c_vtx; 
    for(int i =0; i<nvertices; i++){
        c_vtx=mesh[i];
        if(c_vtx->val == inside_label ){
            x_sum += c_vtx->x;
            y_sum += c_vtx->y;
            z_sum += c_vtx->z;
            n+=1;
        }
    }
    x_sum /= n;
    y_sum /= n;
    z_sum /= n;

    for(int i =0; i<nvertices; i++){
        c_vtx=mesh[i];
        if(c_vtx->val != inside_label ){
            euclidian_distance=sqrt(dist(c_vtx->x, x_sum) + dist(c_vtx->y, y_sum) + dist(c_vtx->z, z_sum));
            c_vtx->dist = c_vtx->dist- euclidian_distance;
        }
    }

return(0);
}


int check_if_border(vertex* c_vtx){
    vertex* t_vtx;

    //update border voxels
    for(int j=0; j<c_vtx->nngh; j++){
        t_vtx=c_vtx->ngh[j];
            if( t_vtx->inside == TRUE){ 
                if(c_vtx->inside != TRUE )
                    //printf("target is inside, current is outside\n");
                    //Current is outside, target is inside, therefore current is border point.
                    return(TRUE);
            }
            else {
                if(c_vtx->inside ==TRUE){
                //Current vertex is inside, target is outside, therefore current is border
                    return(TRUE);
                }   
            }  
    }
    return(FALSE);
}

int add_to_list(void* list, int list_size, int* n, void* element ){
    *n += 1;
    if(list_size==sizeof(double) ){
        double** list_d=list;
        double* element_d= element;
        *list_d=realloc(*list_d, sizeof(**list_d) * *n);
        if( *list_d ==NULL) return(1);
        (*list_d)[*n-1]=*element_d;
    } 

    return(0);
}

double projection(vertex* first_vertex, vertex* second_vertex,vertex* third_vertex, double* x, double* y, double* z){
    //Derivation of shortest distance from first_vertex to edge between second and third vertex:
    //
    //1) D=(x2-x0-t(x1-x0))^2 + (y2-y0-t(y1-y0))^2 + (z2-z0-t(z1-z0))^2
    //2) dD/dt = 2(x2-x0-t(x1-x0))(x1-x0) + 2(y2-y0-t(y1-y0))(y1-y0) + 2(z2-z0-t(z1-z0))(z1-z0) = 0
    //3) 2(x2-x0)(x1-x0)-2*t(x1-x0)(x1-x0) + 2(y2-y0)(y1-y0)-2*t(y1-y0)(y1-y0) + 2(z2-z0)(z1-z0)-2*t(z1-z0)(z1-z0) =0
    //4) (x2-x0)(x1-x0) + (y2-y0)(y1-y0)+ (z2-z0)(z1-z0) =t(x1-x0)(x1-x0) + t(y1-y0)(y1-y0) + t(z1-z0)(z1-z0)
    //
    //      (x2-x0)(x1-x0) + (y2-y0)(y1-y0)+ (z2-z0)(z1-z0)
    //5) t=---------------------------------------------------
    //          (x1-x0)^2 + (y1-y0)^2 + (z1-z0)^2
    //       
    double x0=second_vertex->x;//inside
    double y0=second_vertex->y;
    double z0=second_vertex->z;
    double x1=third_vertex->x;//inside
    double y1=third_vertex->y;
    double z1=third_vertex->z;   
    double x2=first_vertex->x;//outside
    double y2=first_vertex->y;
    double z2=first_vertex->z;   
    double num= (x2-x0)*(x1-x0)+(y2-y0)*(y1-y0)+(z2-z0)*(z1-z0);
    double den= (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0);
    double t=num/den;
    if( t< 0) t=0;
    else if (t >1) t=1; 
    *x=x0+t*(x1-x0);
    *y=y0+t*(y1-y0);
    *z=z0+t*(z1-z0);
    //printf("\tSecond\tThird\tFirst\n");
    //printf("x %f %f, %f --> %f\n", x0, x1, t , *x);
    //printf("y %f %f, %f --> %f\n", y0, y1, t , *y);
    //printf("z %f %f, %f --> %f\n", z0, z1,t , *z);
    return(0);
}

double interpolate_distance(vertex* vertex1, vertex* vertex2, double x, double y, double z){
    double x0=vertex1->x;
    double y0=vertex1->y;
    double z0=vertex1->z;
    double x1=vertex2->x;
    double y1=vertex2->y;
    double z1=vertex2->z; 
    double distance1=vertex1->dist;
    double distance2=vertex2->dist;
    double total_distance=sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
    double part_distance=sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0));
    //when p==p0 -> part_distance==0 --> second_distance==0
    double out=distance1*(1-part_distance/total_distance )+distance2*(part_distance/total_distance);

    return(out);
}

int is_neighbour(vertex* current, vertex* target){
    for(int i=0; i< current->nngh; i++){
       vertex* ngh=current->ngh[i];
       if(ngh == target){
           return(TRUE);
       }
    }
    return(FALSE);
}

double find_smallest_distance(vertex* first_vertex){
    vertex* second_vertex, *third_vertex;
    double* distance_list=NULL; //List of the distances between the 
    double distance;
    double min_distance;
    int ndistances;
    double x, y, z;
    int n_virtual;
    ndistances=0;


    /* 
    #############
    ## Example ##
    ############# 
    v0=first vertex
    [vi] = vertex is "inside" expanding region
    v3=second vertex
    v4=third vertex
  
            v2
          / |  \
       [v3]-v0- v1  
          \ | /
           [v4]

                   v0 
                 _/|\
            9  _/  | \_  8
             _/    |   \_
           _/    5 |     \_
          /        |       \
        [v3] ------------- [v4]
                  v3.5
    [v3] and [v4] are both inside, so we will calculate the distance between v0 projected onto the midpoint
    between [v3] and [v4]. The projected distance, 5, is smaller than the other distances (i.e., 8 and 9),
    so this is taken to be the distance between v0 and the region of inside points.

    */

    //check the neighbouring vertices around the first vertex
    for(int j=0; j < first_vertex->nngh; j++){
        second_vertex=first_vertex->ngh[j];
        //if any of the surroudning vertices are inside
        if(second_vertex->inside == TRUE){
            //printf("Second: %f %f %f\n", second_vertex->x, second_vertex->y, second_vertex->z);
            //then check the neighbours around the second vertex to find another vertex in common with the first one
            for(int i=0; i < second_vertex->nngh; i++){
                third_vertex=second_vertex->ngh[i];
                //If the second and third vertices are both "inside" and are neighbours, calculate the projection
                //from the first vertex to the midpoint on the edge between second and third vertices. 
                if ( third_vertex->inside == TRUE && is_neighbour(first_vertex, third_vertex)==TRUE  ){
                    //Calculate the shortest distance between the first vector and the edge connecting the 
                    //second and third vector
                    projection(first_vertex, second_vertex, third_vertex, &x, &y, &z);
                    //
                    distance=interpolate_distance(second_vertex, third_vertex, x, y, z);
                    
                    distance += sqrt(dist(x, first_vertex->x)+dist(y, first_vertex->y)+dist(z, first_vertex->z));
                    if ( add_to_list(&distance_list, sizeof(distance), &ndistances, &distance ) != 0 ) 
                        pexit("Error: could not allocate memory in <add_to_list>.",  __FILE__, __LINE__, 1 );
                
                    //printf("Projection: %d, %f\n", ndistances, distance );
                }
            }
            
            distance=sqrt(dist(second_vertex->x, first_vertex->x)+dist(second_vertex->y, first_vertex->y)+dist(second_vertex->z, first_vertex->z));

            distance += second_vertex->dist;
            if( add_to_list(&distance_list, sizeof(distance), &ndistances, &distance) != 0 ) 
                pexit("Error: could not allocate memory in <add_to_list>.",  __FILE__, __LINE__, 1 );
            //printf("Edge: %d, %f\n", ndistances, distance );
        }
    }
    
    //Find shortest distance in "distance_list"
    min_distance=distance_list[0];
    for(int i=1; i< ndistances; i++){
        //printf("%d: %f\n", i, distance_list[i]);
        if( min_distance > distance_list[i]) min_distance=distance_list[i];
    }
    //printf("Minimum: %f\n", min_distance);
    free(distance_list);
    return(min_distance);
}


//need to initialize edges

int find_distances(vertex** mesh, int n_vtx,  float dt){
    //This function calculates the distance between vertices based on discrete steps along the edges connecting the vertices together.
    //These steps are referred to as time steps of length dt. The distance of a vertex is the amount of time steps it takes to get to it. 
    float t=0;
    int n_outer_vtx=-1;
    vertex *c_vtx, *t_vtx, *p_vtx;
    edge* e;
    double distance, cp_dist;
    int var;
    int skip_alt=0;
    int n_border_vtx=-1;
    

/*  "|", "/", "_" = edge
    "vi" = node          
    "[vi]" = inside node
        t = 1  
                  _ v3  
                _/  |
              _/    |       
          8 _/     5|        
           /        |
          x         x
        [v1]--------[v2] 

        t = 2  
                  _ v3  
                _/  |
              _/    |       
          8 _/     5|        
           x        x
          |         |
        [v1]--------[v2] 

        .
        .
        .

        t = 5  
                  _[v3]  
                _/  x
              x/    |       
          8 _/     5|        ==> v3 is 5 distance unit away from the "inside" region of [v1] and [v2]
           /        |
          |         |
        [v1]--------[v2] 



    */


    while(n_border_vtx != 0 ){
        printf("%f\t%d\t%d\n",t, n_outer_vtx, n_border_vtx);
        //Set the number of outer vertices to 0
        n_outer_vtx=n_border_vtx=0;
        //increment t by dt (keeps track of total distance travelled along edges)
        t += dt;
        //for all vertiex points
        for(int i=0; i< n_vtx; i++){
            c_vtx=mesh[i]; //current vertex "c_vtx"
            //Count the number of vertices that are not within the region of interest (ROI) (i.e., they are exterior) and
            //are connected by edges to other vertices 
            if (c_vtx->inside != TRUE && c_vtx->n_incoming_edges > 0) n_outer_vtx += 1;
            //Check if c_vtx is an interior, border point
            if( c_vtx->inside == TRUE && c_vtx->border==TRUE){
                var=0;
                //for edges leaving the current vertex
                for(int j=0; j < c_vtx->n_outgoing_edges; j++){
                    e=c_vtx->outgoing_edges[j];
                    t_vtx=e->end; //the target vertex is the one at the end of the current edge
                    //if this edge leads to an outer vertex, check if we have travelled
                    //to this outer vertex
                    if(e->travelled != TRUE && e->end->inside != TRUE ){
                        var += 1;
                        //the total distance to this outer vertex is the
                        //distance travelled to the current vertex plus the edge 
                        //length
                        
                        distance = find_smallest_distance(t_vtx);// c_vtx->dist + e->length;
                        //printf("%f %f\n", distance, t );
                        if (distance <= t){
                            //if t is greater than or equal to distance,
                            //then the edge has been travelled
                            e->travelled=TRUE; //mark edge as travelled
                            if(skip_alt == 1){
                                t_vtx->previous=c_vtx;
                            }
                                //if there are no more edges connecting to the target vertex
                                //then it goes from being outside to inside. 
                                //Its exterior neighbours are also added to the list of frontier points
                            t_vtx->inside=TRUE;
                            t_vtx->border=TRUE;
                            for(int k=1; k < t_vtx->n_incoming_edges; k++){  
                                t_vtx->incoming_edges[k]->travelled=TRUE;
                            } 
                            t_vtx->dist=distance; 
                        }
                    }
                }
                //There are no more untravelled edges for this vertex
                //so it is no longer a border point
                if(var == 0) c_vtx->border=FALSE;
                else n_border_vtx++;
            }
        }
    }
    //for(int i=0; i< n_vtx; i++) if(mesh[i]->dist == MAX ){
    //mesh[i]->dist=0;
    //}
}


/***************************************************************************************************
 *Author: Thomas Funck 
 *Contact: thomas.funck@mail.mcgill.ca
 *Date: September 26, 2014
 *Location: Montreal Neurological Institute, Montreal QC
 *
 *Name: surf_dist  
 *Synopsis: - program to calculate minimum distance from a region to any point on a vertex mesh
 *
 *Inputs: 
 *          1) Mesh composed of vertices in .obj format
 *          2) Binary mask file (.txt) containing "1" for vertices inside region
 *
 *Outputs: 
 *          1) File (.txt) containing distance values at each vertex
 * 
 *
 *Description:
 *  Expands from region at discrete time intervals and calculates edge-wise distance 
 *  from region to each vertex outside the region using projection.
 *  Performance can be improved by providing a higher resolution mesh. 
 *
 *
 *Version History:
 *
 *September 26, 2014
 *   Version-3.0: Compiles and runs correctly with and without option to skip alternate vertices.
 *
 *
 *
 *
 *
 *
 *
 *
 * ****************************************************************/
int main(int argc, char** argv){
    if (argc < 4 || strcmp(argv[0],"-help") ==0  ) useage();
    
    int i=1;
    float dt=0.1;
    if(strcmp(argv[i++],"-dt") == 0){
        dt=atof(argv[i++]);
        printf("%f\n");
    }
    else i=1;

    char* mesh_filename=argv[i++];

    char* vertex_values_filename=argv[i++]; 
    FILE* outputfile=fopen(argv[argc-1], "w+");
    int nvertices, n_border_vtx,n_inner_vtx, n_outer_vtx, n_frontier_vtx=-1;
    int *n_ngh, **ngh;
    int max_vertex_length=4;
    
    char *file_inputs[]={mesh_filename, vertex_values_filename, NULL}; //={mesh_filename, vertex_values_filename};
    float inside_label=1, outside_label=0;
    vertex  **border, **inside, **outside, **frontier;
    VERBOSE=FALSE;
    //Check input files to make sure they exist
    if (check_input_files(file_inputs) != 0) exit(1);
    
    //Read in vertex locations and values
    if(VERBOSE) printf("Reading in vertices.\n");
    vertex** mesh=read_vertices(mesh_filename,  vertex_values_filename, &nvertices, inside_label  );
    if(VERBOSE) printf("Number of vertices: %d\n", nvertices);
    //Find vertex neighbours (ngh)
    if(VERBOSE) printf("Interpolate sphere.\n");
    interpolate_sphere( mesh_filename , &n_ngh, &ngh);
    //Link vertex ngh together
    if(VERBOSE) printf("Linking to neighbours.\n");
    linkVerticesToNeighbours(mesh,  n_ngh, ngh,  nvertices, max_vertex_length);
    find_distances( mesh, nvertices,  dt);


    for(int i=0; i< nvertices; i++) {
        fprintf( outputfile, "%f\n", mesh[i]->dist);
    }
    return 0;
}

void useage(){
    printf("Name:\nsurf_dist ~ calculate distances along mesh surface.\n");
    printf("Description:\nFor a given mask defined on a surface mesh, calculate the minimum distance from thie region to every vertex exterior ");
    printf("to this region on the mesh.\n");
    printf("Useage:\n<-dt> increment> input_mesh.obj vertex_values.txt output.txt\n");
    exit(1);
}
