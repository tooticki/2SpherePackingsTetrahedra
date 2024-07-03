// check conjectured density (except near the conjectured max) for the densest tetrahedron s.t. one sphere is in contact with the three other spheres.
// homothety -> one can assume that two spheres are in contact, say ab
// sliding -> one can assume either a second pair of spheres in contact, say ac or cd (two cases to consider), or a sphere support of radius r
// we thus have four degree of freedom

// Library for interval arithmetic
#include <boost/numeric/interval.hpp>
using namespace boost::numeric;
using namespace interval_lib;

// Threads for parallel programming
#include <thread>         // std::thread
int threads_finished = 0;


#include <array>
#include <iostream>


// Fix rounding policies for the transcendental functions
typedef interval<double, policies<save_state<rounded_transc_std<double>>, checking_base<double>>> I;

typedef std::array<I,6> block; // 5 pour gagner 1/6 de mémoire ? 96 bits, i.e. 8bits for each I
 // maximum nb of blocks in memory
//const int maxactifs=20000000; //32% mem
const int maxactifs=50000000;

// size ratio for halite
I r=sqrt(I(2))-I(1);
I u=I(1);

// radii
I ra=u;
I rb=u;
I rc=u;
I rd=u;



/******************/
/* case selection */
/******************/

// radii used in compilation to optimize the formula for the radius of sphere support
#define uuuu
// TODO : lister les cas

// contact_ac or contact_cd for a second contact along ac or cd
// for a support sphere of size r
#define support_sphere_r

/******************/


void print_block(block B)
{
    printf("(");
    for(int i=0;i<6;i++)
        printf("\nRIF(%.50f,%.50f),",lower(B[i]),upper(B[i]));
    printf("\b \n)\n");
}

// radius of the sphere support of a tetrahedron
#include "radius_general.cpp"
// length of ac as a function of the other edge length and a radius r
#include "ac.cpp"
// other routines for tetrahedra (volume, angles, density etc)
#include "routines.cpp"

// lower bound on the conjectured maximal density
#if defined (uuuu)
    const I delta_bound=I(779635570044252)/I(1000000000000000);
//const I delta_bound=I(784688454045207)/I(1000000000000000); // rrrr with two incident edges rr stretched
#else //TODO
    const I delta_bound=I(812542027810834)/I(1000000000000000); // r111 with one edge 11 stretched
#endif
// stretched edges
const I t=I(4)*sqrt(I(2)*r/(I(1)+I(2)*r)); // 11 in 111r
const I s=r*sqrt(I(2)*sqrt(I(6))+I(6)); // rr in rrrr

// track densest deleted block
block densest_block;
float delta_max=0;





/*************************************************************/
/*************************************************************/
/*************************************************************/

// decide whether to keep block B for further refinement)
bool keep(block B)
{
    #if defined(support_sphere_r)
    B[1]=ac(B);
    if (empty(B[1])) return false; // not such FM-tetrahedra
    #endif

    //symmetries
//    #if defined(ruuu) // one can assume CD>=BC>=BD
//    if (upper(B[5])<=lower(B[3]) || upper(B[3])<lower(B[4])) return false;
//    #endif

    // skip block if it is within epsilon from the conjectured densest tetrahedron
    // (epsilon and the edge lengths depend on the considered case)
    #if defined(uuuu)
    I eps=I(1)/I(80);
    int i=0;
    for(int j=0;j<6;j++) // count contacts
        if (upper(B[j])<lower(I(2)+eps)) i+=1;
    if (i==6) return false; // 6 eps-contacts -> skip
    #elif defined(ruuu)
    I eps=I(1)/I(499); // prouvé
//    I eps=I(1)/I(80); // pas prouvé (essai pour voir si ça passe)
    int i=0;
    for(int j=0;j<6;j++) // count contacts
        if (upper(B[j])<lower(I(2)+eps)) i+=1;
    // 5 contacts and one stretched edge 11 -> skip
    if (i==5 && (upper(abs(B[3]-t))<lower(eps) || upper(abs(B[4]-t))<lower(eps) || upper(abs(B[5]-t))<lower(eps))) return false;
    #endif

    // reject blocks which do not correspond to any tetrahedron
    if (!tetrahedral(B[0],B[1],B[2],B[3],B[4],B[5]))
        return false;

    // radius of the support sphere of the block B
    #if defined(support_sphere_r)
        I rhoB=r;
    #else
    I rhoB=radius(B[0],B[1],B[2],B[3],B[4],B[5]);
    if (lower(rhoB)>upper(r)) return false; // too large -> reject block
    #endif
    
    // density of the block B
    I deltaB=density(B[0],B[1],B[2],B[3],B[4],B[5]);

    // if the block has density surely lower than the bound, reject it
    if (upper(deltaB)<lower(delta_bound))
    {
        // update the lower bound on the density of the densest rejected block
        // actually we should check that the support sphere has radius at most r...
        if (lower(deltaB)>delta_max)
        {
            densest_block=B;
            delta_max=lower(deltaB);
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
#if defined(support_sphere_r)
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

void bound_density_in_block(block B){
    printf("thread id:%i\n", std::this_thread::get_id());
    printf(" Bound density in block: ");
    print_block(B);
    // blocks stores all the active blocks (only one at the beginning)
    block* blocks=new block[1];
    blocks[0]=B;
    int actifs=1; // counts the active blocks
    int newactifs=0; // counts the active blocks created by halving during one step
    int del=0;
    int step=0;

    while (actifs>0) // split blocks step by step while there are still some
    {
        int newdel=0;
        if (actifs>=0) printf(" step %2d: %9d blocks considered",step,2*actifs); fflush(stdout);
        // to store the active blocks created by halving during this step
        block* newblocks=new(std::nothrow) block[2*actifs]; // each block will give at most two new blocks
        newactifs=0;
        int timer=1 ;

        for(int i=0;i<actifs;i++)
        {
	    // if (actifs>=10 && i*10/actifs==timer) {timer++;std::cout << "." << std::flush;}

            // BB is the pair of blocks obtained by halving blocks[i]
            auto BB=split(blocks[i]);

            // we put in newblocks each block in the pair BB which has to be kept
            if (keep(BB.first))  {newblocks[newactifs]=BB.first;  newactifs++;} else {newdel++;}
            if (keep(BB.second)) {newblocks[newactifs]=BB.second; newactifs++;} else {newdel++;}
        }
        del+=newdel;
        if (actifs>=0) printf("|%9d| blocks deleted (%2d\%)\n",newdel,100*newdel/2/actifs);

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
/*            I V=vol(B[0],B[1],B[2],B[3],B[4],B[5]);
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
            printf("D=RIF(%f,%f)\n",lower(z),upper(z)); */
            }
            break;
        }
        step++;
    }
    delete blocks;

    printf("Total number of considered blocks: %d\n",del);

    printf("Block with the highest lower bound on the density (%.20f):\n",delta_max);
    print_block(densest_block);
    threads_finished+=1;
    printf("%i threads terminated", threads_finished);
}


int main(int argc, char *argv[])
{
    using namespace std;
    printf("Tetrahedron ");
    if (ra==u) printf("1"); else printf("r");
    if (rb==u) printf("1"); else printf("r");
    if (rc==u) printf("1"); else printf("r");
    if (rd==u) printf("1"); else printf("r");
    printf("\n");

    I gap=I(2)*hull(I(0),r);


    // initial block
    block B;
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

    printf("Initial block: ");
    print_block(B);


    int sn = 6;
    int thread_number = pow(2,sn);
    // Subdivide the initial block into 2^sn
    block* initial_blocks=new(std::nothrow) block[thread_number]; // each block will give at most two new blocks	    
    initial_blocks[0] = B;
    block* new_blocks=new(std::nothrow) block[thread_number]; // each block will give at most two new blocks	    
    vector<thread> threads(thread_number);
    
    // generate 2^sn sub blocks by halving 4 times
    for (int i = 0; i<sn; i++){
	printf("%i",i);
	//new_blocks = {};
	for(int j = 0; j<pow(2,i); j++){
	    printf(" %i",j);
	    auto bb = split(initial_blocks[j]);
	    new_blocks[2*j] = bb.first;
	    new_blocks[2*j+1] = bb.second;	    
	}
	for(int j = 0; j<pow(2,i+1); j++){
	    printf(" %i",j);
	    initial_blocks[j] = new_blocks[j];
	}
    }

    // spawn 16 threads on 16 blocks:
    for (int i = 0; i < thread_number; i++) {
	print_block(initial_blocks[i]);
        threads[i] = thread(bound_density_in_block, initial_blocks[i]);
	printf("%i",i);
    }

    for (auto& th : threads) {
        th.join();
    }
    

    return 0;
}
