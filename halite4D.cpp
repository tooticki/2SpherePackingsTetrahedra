// check conjectured density (except near the conjectured max) for the densest tetrahedron s.t. one sphere is in contact with the three other spheres.
// homothety -> one can assume that two spheres are in contact, say ab
// sliding -> one can assume either a second pair of spheres in contact, say ac or cd (two cases to consider), or a sphere support of radius r
// we thus have four degree of freedom

#include <iostream>
using namespace std;

// Library for interval arithmetic
#include <boost/numeric/interval.hpp>
using namespace boost::numeric;
using namespace interval_lib;

// Boost asio for parallel programming
#include <boost/bind.hpp>
#include <boost/thread.hpp>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>

// Timer
#include <chrono>
using namespace std::chrono;

// Fix rounding policies for the transcendental functions
typedef interval<double, policies<save_state<rounded_transc_std<double>>, checking_base<double>>> I;

typedef std::array<I,6> block; // 5 pour gagner 1/6 de mémoire ? 96 bits, i.e. 8bits for each I

// Maximum nb of blocks in memory
const int maxactifs=10000000;

// Sphere radii for halite
I r=sqrt(I(2))-I(1);
I u=I(1);

/****************************** Case selection ******************************/

// Radii used in compilation to optimize the formula for the radius of sphere support
#define uuuu
// TODO : lister les cas

// `contact_ac` or `contact_cd` for a second contact along ac or cd
// `support_sphere_r` for a support sphere of size r
#define support_sphere_r


// Only consider tetrahedra close to the degenerate ones from uuuu support_sphere_r: the one with at least one stretched edge
// in eps_suuuu - neighborhood
#define only_degenerate

/************************************************************/

// Radii values of the tetrahedra and
// lower bounds on the conjectured maximal density
#if defined (uuuu)
I ra=u, rb=u, rc=u, rd=u;
#elif defined (uuur)
I ra=u, rb=u, rc=u, rd=r;
#elif defined (uurr)
I ra=u, rb=u, rc=r, rd=r;
#elif defined (urur)
I ra=u, rb=r, rc=u, rd=r;
#elif defined (rruu)
I ra=r, rb=r, rc=u, rd=u;
#elif defined (urrr)
I ra=u, rb=r, rc=r, rd=r;
#elif defined (rrru)
I ra=r, rb=r, rc=r, rd=u;
#elif defined (rrrr)
I ra=r, rb=r, rc=r, rd=r;
#elif defined (ruuu)
I ra= r, rb=u, rc=u, rd=u;
#else
I ra = I::empty(), rb=I::empty(), rc=I::empty(), rd=I::empty();
#endif

// Lower bounds on the conjectured density and the conjectured optimal blocks
#if defined (uuuu)
const I delta_bound=I(779635570044252)/I(1000000000000000); // 1111 tight
#elif defined (uurr) || defined(urur) || defined(rruu)
const I delta_bound=I(810466032832064)/I(1000000000000000); // 11rr tight
#elif defined (urrr) || defined(rrru)
const I delta_bound=I(806503318194798)/I(1000000000000000); // 1rrr tight
#elif defined(ruuu) ||  defined (uuur) 
const I stretched_edge_length=I(4)*sqrt(I(2)*r/(I(1)+I(2)*r)); //  stretched edge 11 in 111r
const I delta_bound=I(812542027810834)/I(1000000000000000); // r111 with one edge 11 stretched
#elif defined (rrrr)
const I stretched_edges_length=r*sqrt(I(2)*sqrt(I(6))+I(6)); //  stretched edges rr in rrrr
const I delta_bound=I(784688454045207)/I(1000000000000000); // rrrr with two incident edges rr stretched
#else const I delta_bound = I::empty();
#endif

block B_tight = {ra+rb, ra+rc, ra+rd, rb+rc, rb+rd, rc+rd};

// Values of eps obtained with the local maximality proof 
#if defined (uuuu)
I eps_loc = I(1)/I(80);
#elif defined (uuur) || defined(ruuu)
I eps_loc = I(1)/I(499);
#elif defined (uurr) || defined(urur) || defined(rruu)
I eps_loc = I(1)/I(317);
#elif defined (urrr) || defined(rrru)
I eps_loc = I(1)/I(410);
#elif defined (rrrr)
I eps_loc = I(1)/I(243);
#else
I eps_loc = I::empty();
#endif

// Value of eps needed to exclude the degenerate case when ac cannot be computed 
#if defined(uuuu) && defined (support_sphere_r)
I eps_suuuu=I(3)/I(10); //1/17 4s, 1/18 ??? 
#endif

// Functions to verify if the block is in eps-neighborhood of the optimum
#if defined (uuuu) ||  defined (uurr) || defined(urur) || defined(rruu) ||  defined (urrr) || defined(rrru)
bool is_eps_optimal(block B){
    int contacts=0;
    for(int i=0;i<6;i++) // count contacts
        if (upper(B[i]) < lower(B_tight[i]+eps_loc)) contacts++;
    return contacts==6;
}
#elif  defined (uuur) || defined(ruuu)
bool is_eps_optimal(block B){
    int contacts=0, stretched=0;
    for(int i=0;i<6;j++) {// count contacts and stretched edges
        if (upper(B[i]) < lower(B_tight[i] + eps_loc)) contacts++;
#if defined(uuur)
	if ((i==1||i==3) && lower(B[i]) > upper(stretched_edge_length - eps_loc))   stretched++;
#elif defined(ruuu)
	if ((i>=3) && lower(B[i]) > upper(stretched_edge_length - eps_loc))   stretched++;
#endif
    }    
    return contacts==5 && stretched==1;
}
#elif defined (rrrr)
bool is_eps_optimal(block B){
    int contacts=0, stretched=0;
    for(int i=0;i<6;j++) {// count contacts and stretched edges
        if (upper(B[i]) < lower(B_tight[i] + eps_loc)) contacts++;
	if (lower(B[i]) > upper(stretched_edges_length - eps_loc)) stretched++;
    }    
    return contacts==4 && stretched==2;
}
#endif

// radius of the support sphere of a tetrahedron
#include "radius_optimized.cpp"
// length of ac as a function of the other edge length and a radius r optimized for radii
#include "ac_optimized.cpp"
// other routines for tetrahedra (volume, angles, density etc)
#include "routines.cpp"

block* initial_blocks;

// Track densest deleted block and thel lower bound of its density
block densest_block;
float delta_max=0;
// The mutex to synchronise `denstest_block` and `delta_max` which are shared between threads
boost::mutex densest_blockMutex;

/************************************************************/
/************************************************************/
/************************************************************/

void print_block(block B){
    printf("(");
    for(int i=0;i<6;i++)
        printf("\nRIF(%.50f,%.50f),",lower(B[i]),upper(B[i]));
    printf("\b \n)\n");
}

// decide whether to keep block B for further refinement)
bool keep(block B){
    //print_block(B);
#if defined(support_sphere_r) && !defined(only_degenerate)
    B[1]=ac(B, false); // compute edge ac
    if (empty(B[1]))  return false;  // no such FM-tetrahedra
#endif
    // ********** symmetries
    //    #if defined(ruuu) // one can assume CD>=BC>=BD
    //    if (upper(B[5])<=lower(B[3]) || upper(B[3])<lower(B[4])) {print_block(B); return false;} 
    //    #endif

    // skip block if it is within epsilon from the conjectured optimum
    // (epsilon and the edge lengths depend on the considered case)
    if (is_eps_optimal(B))   return false;

    // skip particular cases treated separately: when support sphere is of radius r and can not compute ac
#if defined(uuuu) && defined(support_sphere_r)  && !defined(only_degenerate)
    int contacts = 0, stretched = 0;
    for(int j=0;j<6;j++) {// count contacts
        if (upper(B[j])<lower(I(2)+eps_suuuu)) contacts++;
	if (lower(B[j])>upper(I(2)+I(2)*r-eps_suuuu)) stretched++;  }
    if (contacts >= 4 && stretched >= 1)  return false;  // 4 eps-tight and 1 eps-stretched edges
#endif
    
    // reject blocks which do not correspond to any tetrahedron
    if (!tetrahedral(B[0],B[1],B[2],B[3],B[4],B[5])) return false;

    // radius of the support sphere of the block B
#if defined(support_sphere_r)
    I rhoB=r;
#else
    I rhoB=radius(B[0],B[1],B[2],B[3],B[4],B[5]);
    if (lower(rhoB)>upper(r))  return false;  // support sphere is too large -> reject block
#endif
    
    // density of the block B
    I deltaB=density(B[0],B[1],B[2],B[3],B[4],B[5]);
    //cout << lower(deltaB)<<endl;
    
    // if the block has density surely lower than the bound, reject it
    if (upper(deltaB)<lower(delta_bound)) {
	// update the lower bound on the density of the densest rejected block
	// ! remove it if too slow !
	if (lower(deltaB)>delta_max) {
	    densest_blockMutex.lock();
	    densest_block=B;
	    delta_max=lower(deltaB);
	    densest_blockMutex.unlock();
	}
	return false;
    }
    
    // now the density MAY be larger than the conjectured upper bound
    // if the radius of the support sphere is not SURELY at most r then the block contains both FM- and non-FM-tetrahedra: -> we need to refine
#if !defined(support_sphere_r)
    if (upper(rhoB)>=lower(r)) return true;
#endif

    // from now on the block as SURELY radius at most r, i.e. it contains only FM-tetrahedra
    // we can consider density: it is surely more than the bound it yields a counter example!
    if (lower(deltaB)>upper(delta_bound))
    {
        printf("The following block has density at least %f:\n",lower(deltaB));
        print_block(B);
        printf("r=%f,%f\n",lower(rhoB),upper(rhoB));
        I V=vol(B[0],B[1],B[2],B[3],B[4],B[5]);
        I C=cov(B[0],B[1],B[2],B[3],B[4],B[5]);
        printf("vol=%f,%f\n",lower(V),upper(V));
        printf("cov=%f,%f\n",lower(C),upper(C));
        throw "denser block found!";
    }

    // the block contains only FM-tetrahedra BUT its density interval intersects delta_bound -> we need to refine
    return true;

}

// split the block B along its largest edge
std::pair<block,block> split(block B)
{
    int i=0;
    for(int j=1;j<6;j++){
	//ifi the sphere support has radius r, ac is derived from other edges so ignore it
#if defined(support_sphere_r) && !defined(only_degenerate)
	if (j==1) continue;
#endif
        if (width(B[j])>width(B[i])) i=j;}
    std::pair<I,I> e=bisect(B[i]);
    block C=B;
    B[i]=e.first;
    C[i]=e.second;
    return (std::pair<block,block>(B,C));
}


/*************************************************************/
/*************************************************************/
/*************************************************************/

int bound_density_in_block(int i, bool verbose=false){
    block B = initial_blocks[i];
    // blocks stores all the active subblocks (only one at the beginning)
    block* blocks=new block[1];
    blocks[0]=B;
    int actifs=1; // counts the active blocks
    int newactifs=0; // counts the active blocks created by halving during one step
    int del=0;
    int step=0;

    while (actifs>0) // split blocks step by step while there are still some
	{
	    int newdel=0;
	    if (verbose && actifs>=0) printf(" step %2d: %9d blocks considered",step,2*actifs); fflush(stdout);
	    // to store the active blocks created by halving during this step
	    block* newblocks=new(std::nothrow) block[2*actifs]; // each block will give at most two new blocks
	    newactifs=0;
	    int timer=1 ;

        for(int i=0;i<actifs;i++)
        {
	     if (verbose && actifs>=10 && i*10/actifs==timer) {timer++;std::cout << "." << std::flush;}

            // BB is the pair of blocks obtained by halving blocks[i]
            auto BB=split(blocks[i]);

            // we put in newblocks each block in the pair BB which has to be kept
            if (keep(BB.first))  {newblocks[newactifs]=BB.first;  newactifs++;} else {newdel++;}
            if (keep(BB.second)) {newblocks[newactifs]=BB.second; newactifs++;} else {newdel++;}
        }
        del+=newdel;
        if (verbose && actifs>=0) printf("|%9d| blocks deleted (%2d\%)\n",newdel,100*newdel/2/actifs);
        // newblocks will play the role of blocks in the next step
        delete blocks;
        blocks=newblocks;
        actifs=newactifs;        
        
        if (actifs>=maxactifs)
        {
            printf("To much blocks! Give up:(\n");
            printf("Random remaining blocks:\n");
            for(int j=0;j<3;j++){
            block B=blocks[rand()%actifs];
            print_block(B);
            I V=vol(B[0],B[1],B[2],B[3],B[4],B[5]);
            I C=cov(B[0],B[1],B[2],B[3],B[4],B[5]);
            I d=density(B[0],B[1],B[2],B[3],B[4],B[5]);
            printf("vol=%f,%f\n",lower(V),upper(V));
            printf("cov=%f,%f\n",lower(C),upper(C));
            printf("density=%f,%f\n",lower(d),upper(d));
            I ab=B[0];
            I ac=B[1];
            I ad=B[2];
            I bc=B[3];
            I bd=B[4];
            I cd=B[5];
            I z=solid(ab,ac,ad,bc,bd,cd);
            printf("A=RIF(%f,%f)\n",lower(z),upper(z));
            z=solid(bc,bd,ab,cd,ac,ad);
            printf("B=RIF(%f,%f)\n",lower(z),upper(z));
            z=solid(cd,ac,bc,ad,bd,ab);
            printf("C=RIF(%f,%f)\n",lower(z),upper(z));
            z=solid(ad,bd,cd,ab,ac,bc);
            printf("D=RIF(%f,%f)\n",lower(z),upper(z)); 
            }
            break;
        }
        step++;
    }
    delete blocks;
    return del;
    //  printf("Total number of considered blocks: %d\n",del);
}

void split_and_run_block(block B, int N=10){
    printf("Initial block: ");
    print_block(B);
    
    // Parallelization: we subdivide the initial block into N^4 small blocks (each free edge is divided by N),
    //  we create M threads, and we give the blocks to threads so that each thread is occupied at each moment
#if defined(contact_ac) || defined(contact_cd) || (defined(support_sphere_r) && !defined(only_degenerate)) // only 3 free variables
    const int blocks_number = pow(N,4);
#else
    const int blocks_number = pow(N,5);
#endif
    //int threads_number = 20;
    initial_blocks=new(std::nothrow) block[blocks_number];
    
    //TODO make a version for the case ac (only 3 edges) and support sphere
    float e1_l = lower(B[1]),  e1_r = upper(B[1]);
    float e2_l = lower(B[2]),  e2_r = upper(B[2]);
    float e3_l = lower(B[3]),  e3_r = upper(B[3]);
    float e4_l = lower(B[4]),  e4_r = upper(B[4]);
    float e5_l = lower(B[5]),  e5_r = upper(B[5]);
    I e0 = B[0], e1, e2, e3, e4,e5;
    int block_index = 0;
    printf("Dividing initial block into %i blocks...", blocks_number);
    
    for (int i1=0; i1<N; i1++){	
#if defined(contact_ac) || (defined(support_sphere_r) &&  !defined(only_degenerate) )
	e1 = hull(e5_l + (e5_r-e5_l)*i1/N, e5_l+(e5_r-e5_l)*(i1+1)/N);
#else
	e1 = hull(e1_l + (e1_r-e1_l)*i1/N, e1_l+(e1_r-e1_l)*(i1+1)/N);
#endif
	for (int i2=0; i2<N; i2++){
	    e2 = hull(e2_l + (e2_r-e2_l)*i2/N, e2_l+(e2_r-e2_l)*(i2+1)/N);
	    for (int i3=0; i3<N; i3++){
		e3 = hull(e3_l + (e3_r-e3_l)*i3/N, e3_l+(e3_r-e3_l)*(i3+1)/N);
		for (int i4=0; i4<N; i4++){
		    e4 = hull(e4_l + (e4_r-e4_l)*i4/N, e4_l+(e4_r-e4_l)*(i4+1)/N);
#if defined(contact_cd)
		    initial_blocks[block_index] = {e0, e1, e2, e3, e4, B[5]};
		    block_index++;
#elif defined(contact_ac) || (defined(support_sphere_r) && !defined(only_degenerate))		    
		    initial_blocks[block_index] = {e0, B[1], e2, e3, e4, e1};
		    block_index++;
#else
		    for (int i5=0; i5<N; i5++){
			e5 = hull(e5_l + (e5_r-e5_l)*i5/N, e5_l+(e5_r-e5_l)*(i5+1)/N);
			initial_blocks[block_index] = {e0, e1, e2, e3, e4, e5};
			block_index++;
		    }
#endif
		}			    
	    }
	}
    }

    printf(" Done \n");
    // async leaks memory: so can not divide more than by 10 each edge.... 
    // run_tasks([](int i) { return bound_density_in_block(i); }, 20, blocks_number);
    
    // another solution: boost thread pool
    boost::asio::thread_pool th_pool(20);
    
   
    auto start = high_resolution_clock::now();
    bool verbose = N<=2;
    for (size_t i = 0; i < blocks_number; i++) {
	post(th_pool, boost::bind(bound_density_in_block,  i, false));//[j=i](){cout << "bind"; return boost::bind(bound_density_in_block,  j);});
    }
    th_pool.join();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    printf("Block with the highest lower bound on the density (%.20f):\n",delta_max);
    print_block(densest_block);
    cout << "Number of blocks : "	 << blocks_number  << endl;    
    cout << "Time taken : "   << duration.count()/1000000 << " seconds" << endl << endl;    
}

int main(int argc, char *argv[])
{
    using namespace std;
    using namespace boost;
    using namespace asio;
    printf("Tetrahedron ");
    if (ra==u) printf("1"); else printf("r");
    if (rb==u) printf("1"); else printf("r");
    if (rc==u) printf("1"); else printf("r");
    if (rd==u) printf("1"); else printf("r");
    printf("\n");

    I gap=I(2)*hull(I(0),r);    
    block B;
    
#if defined(support_sphere_r) && defined(uuuu) && defined(only_degenerate) // If we only consider degenerate tetrahedra: also subdivide ac
    cout << "Only degenerate" << endl;
    I eps_gap = hull(I(0), eps_suuuu);
    B[0]=ra+rb; // always contact by homothety
    B[1]=ra+rc+gap;
    B[2]=ra+rd+eps_gap;
    B[3]=rb+rc+eps_gap;
    B[4]=rb+rd+eps_gap;
    B[5]=rc+rd+eps_gap;
    for(int i = 2; i<=5; i++){ // each of the last 4 edges can be stretched and there is no symmetry
	I e = B[i];
	B[i] += I(2)*r - eps_gap;
	split_and_run_block(B,20);
	B[i] = e;
    }
#else
    B[0]=ra+rb; // always contact by homothety
    B[1]=ra+rc+gap;
    B[2]=ra+rd+gap;
    B[3]=rb+rc+gap;
    B[4]=rb+rd+gap;
    B[5]=rc+rd+gap;
    // dimension reduction by sliding
#if defined(contact_ac) // a second contact along ac
    B[1]=ra+rc;
#elif defined(contact_cd) // or along cd
    B[5]=rc+rd;
#elif defined(support_sphere_r) // or only one contact but a support sphere of radius r
    // ac va devoir être recalculée en fonction des autres longueurs d'arêtes
#endif
    split_and_run_block(B,10);
#endif
    return 0;
}
