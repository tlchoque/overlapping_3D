
void segmentedSphere(Delaunay & dt,double xc,double yc, double zc, double r){
	for(Finite_cells_iterator ci = dt.finite_cells_begin(); ci != dt.finite_cells_end(); ci++){
		double r2 =r*r;
		dt.tetrahedron(ci)[0];
		if( ( pow(dt.tetrahedron(ci)[0].x()-xc, 2) + pow(dt.tetrahedron(ci)[0].y()-yc , 2 ) + pow(dt.tetrahedron(ci)[0].z()-zc , 2 ) < r2 ) && 
			( pow(dt.tetrahedron(ci)[1].x()-xc, 2) + pow(dt.tetrahedron(ci)[1].y()-yc , 2 ) + pow(dt.tetrahedron(ci)[1].z()-zc , 2 ) < r2 ) && 
			( pow(dt.tetrahedron(ci)[2].x()-xc, 2) + pow(dt.tetrahedron(ci)[2].y()-yc , 2 ) + pow(dt.tetrahedron(ci)[2].z()-zc , 2 ) < r2 ) &&
			( pow(dt.tetrahedron(ci)[3].x()-xc, 2) + pow(dt.tetrahedron(ci)[3].y()-yc , 2 ) + pow(dt.tetrahedron(ci)[3].z()-zc , 2 ) < r2 )){
			ci->info().m_original_label = 1;
			ci->info().m_label = 1;
		}	
		else {
			ci->info().m_original_label = 0;
			ci->info().m_label = 0;
		}
	}
}

void segmentedCube(Delaunay & dt, double x, double y,double z,double side){
	double left,right,up,down,back,forward;
	left=x;right=x+side;up=y;down=y-side;back=z;forward=z+side;
	for(Finite_cells_iterator ci = dt.finite_cells_begin(); ci != dt.finite_cells_end(); ci++){
		if( (dt.tetrahedron(ci)[0].x() < right && dt.tetrahedron(ci)[0].x() > left && dt.tetrahedron(ci)[0].y() < up && dt.tetrahedron(ci)[0].y() > down && 
			 dt.tetrahedron(ci)[0].z() < forward && dt.tetrahedron(ci)[0].z() > back) &&
			(dt.tetrahedron(ci)[1].x() < right && dt.tetrahedron(ci)[1].x() > left && dt.tetrahedron(ci)[1].y() < up && dt.tetrahedron(ci)[1].y() > down && 
			 dt.tetrahedron(ci)[1].z() < forward && dt.tetrahedron(ci)[1].z() > back) &&
			(dt.tetrahedron(ci)[2].x() < right && dt.tetrahedron(ci)[2].x() > left && dt.tetrahedron(ci)[2].y() < up && dt.tetrahedron(ci)[2].y() > down && 
			 dt.tetrahedron(ci)[2].z() < forward && dt.tetrahedron(ci)[2].z() > back) &&
			(dt.tetrahedron(ci)[3].x() < right && dt.tetrahedron(ci)[3].x() > left && dt.tetrahedron(ci)[3].y() < up && dt.tetrahedron(ci)[3].y() > down && 
			 dt.tetrahedron(ci)[3].z() < forward && dt.tetrahedron(ci)[3].z() > back)){
			ci->info().m_label = 2;
			ci->info().m_original_label = 2;
		}
		//else ci->info().m_label = 0;//test manifold
	}
}

void segmentedPlane(Delaunay & dt, double b, double c,double d){
	for(Finite_cells_iterator ci = dt.finite_cells_begin(); ci != dt.finite_cells_end(); ci++){
		if( ( dt.tetrahedron(ci)[0].y() - b*dt.tetrahedron(ci)[0].x() - c*dt.tetrahedron(ci)[0].z() - d >0 ) &&
			( dt.tetrahedron(ci)[1].y() - b*dt.tetrahedron(ci)[1].x() - c*dt.tetrahedron(ci)[1].z() - d >0 ) &&
			( dt.tetrahedron(ci)[2].y() - b*dt.tetrahedron(ci)[2].x() - c*dt.tetrahedron(ci)[2].z() - d >0 ) &&
			( dt.tetrahedron(ci)[3].y() - b*dt.tetrahedron(ci)[3].x() - c*dt.tetrahedron(ci)[3].z() - d >0 ) ){
			ci->info().m_label = 3;
			ci->info().m_original_label = 3;
		}
	}
}

void segmentation(Delaunay & dt, int windowWidth, int windowHeight, int windowDepth){
	double xc,yc,zc,r;
	xc=(double)( rand()%windowWidth );
	yc=(double)( rand()%windowHeight );
	zc=(double)( rand()%windowDepth );
	r=(double)( rand()%(int)floor(windowHeight/3)   ) + (windowHeight/5);

	xc=650;//650
	yc=450;//450
	zc=250;
	r=250;//250
	segmentedSphere(dt,xc,yc,zc,r);

	double x,y,z,side; // left upper deeper point of square (x,y,z) and volume of cube = side^3
	side=(double)( rand()%(int)floor(windowHeight/4) ) + (windowHeight/3);
	x=(double)(rand()%windowWidth - (int)(side/2) );
	y=(double)(rand()%windowHeight + (int)(side/2) );
	z=(double)(rand()%windowDepth + (int)(side/2) );

	x=0;
	y=450;
	z=0;
	side=450;
	//test non manifold
	/*x=200;
	y=650;
	z=200;
	side=450;*/
	segmentedCube(dt,x,y,z,side);

	double b,c,d;
	b=((double) rand() / (RAND_MAX));
	c=((double) rand() / (RAND_MAX));
	d=windowHeight - (double)( rand()%(windowHeight/4) ) -  windowHeight/4;
	b=0.3;
	c=0.4;
	d=500;
	segmentedPlane(dt,b,c,d);
}