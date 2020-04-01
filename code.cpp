EdgeRecord::EdgeRecord(EdgeIter& _edge) : edge(_edge) {
  // TODO: (meshEdit)
  // Compute the combined quadric from the edge endpoints.
  // -> Build the 3x3 linear system whose solution minimizes the quadric error
  //    associated with these two endpoints.
	Matrix4x4 quadricE = edge->halfedge()->vertex()->quadric + edge->halfedge()->twin()->vertex()->quadric;
	Matrix3x3 A;
	for (size_t x = 0; x <= 2; x++)
	{
		for (size_t y = 0; y <= 2; y++)
		{
			A(x, y) = quadricE(x, y);
		}
	}
	Vector3D w(quadricE(0, 3), quadricE(1, 3), quadricE(2, 3));
	Vector3D b = -w;
  // -> Use this system to solve for the optimal position, and store it in
  //    EdgeRecord::optimalPoint.
	if (abs(A.det()) < 0.000001) {
		optimalPoint = edge->centroid();
	}
	else
	{
		optimalPoint = A.inv()*b;
	}
  // -> Also store the cost associated with collapsing this edg in
  //    EdgeRecord::Cost.
	Vector4D u = Vector4D(optimalPoint, 1);
	score = dot((quadricE*u), u);
}

void MeshResampler::downsample(HalfedgeMesh& mesh) {
	for (auto f = mesh.facesBegin(); f != mesh.facesEnd(); f++)
	{
		double d = -dot(f->normal(), f->halfedge()->vertex()->position);
		Vector4D v = Vector4D(f->normal(), d);
		f->quadric = outer(v, v);
	}
 
	for (auto v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++)
	{
		auto adjFs = v->AdjFaces();
		v->quadric.zero();
		for (auto f : adjFs)
			v->quadric += f->quadric;
	}
 
	MutablePriorityQueue<EdgeRecord> queue;
	for (auto e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++)
	{
		e->record = EdgeRecord(e);
		queue.insert(e->record);
	}
 
	size_t targetNum = mesh.nFaces() / 4;
	while (mesh.nFaces()>targetNum)
	{
		EdgeRecord eR = queue.top();
		queue.pop();
 
		{//remove adjEs' record in queue
			auto adjEs = eR.edge->AdjEdges();
			for (auto adjE : adjEs)
				queue.remove(adjE->record);
		}
 
		VertexIter newV = mesh.collapseEdge(eR.edge);
		if (!mesh.IsValid(newV, "downsample : collapse an edge fail"))
			return;
 
		newV->position = eR.optimalPoint;
 
		{// set adjFs' and newV's quadric
			newV->quadric.zero();
			auto adjFs = newV->AdjFaces();
			for (auto adjF : adjFs) {
				double d = -dot(adjF->normal(), adjF->halfedge()->vertex()->position);
				Vector4D v = Vector4D(adjF->normal(), d);
				adjF->quadric = outer(v, v);
				newV->quadric += adjF->quadric;
			}
		}
 
		{// set adjVs' quadric
			auto adjVs = newV->AdjVertices();
			for (auto adjV : adjVs) {
				adjV->quadric.zero();
				auto adjFs = adjV->AdjFaces();
				for (auto adjF : adjFs)
					adjV->quadric += adjF->quadric;
			}
		}
 
		{// set adjEs' record
			auto adjEs = newV->AdjEdges();
			for (auto adjE : adjEs) {
				adjE->record = EdgeRecord(adjE);
				queue.insert(adjE->record);
			}
		}
	}
}