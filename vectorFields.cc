#include "vectorFields.hh"
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"

#include <string>
#include <queue>
#include <algorithm>

// Helper functions:
ACG::Matrix3x3d orthogonalize(ACG::Vec3d u, ACG::Vec3d v ) {
    ACG::Vec3d Rx = u.normalize();

    ACG::Vec3d Ry = v;
    Ry = Ry - (Ry | Rx) * Rx;
    Ry = Ry.normalize();
 
    ACG::Vec3d Rz = Rx % Ry;
    Rz = Rz.normalize();

    ACG::Matrix3x3d R;
    R = R.fromColumns(Rx, Ry, Rz);

    // std::cout << u << " " << v << " " << R << std::endl;

    return R;
}

ACG::Vec3d VectorFields::transport(ACG::Vec3d w0, OpenMesh::SmartFaceHandle& fh_i, OpenMesh::SmartFaceHandle& fh_j) {

    // Grab vertices a, b, c, and d of the two adjacent
    // triangles i and j according to the following labels:
    //                
    //                         b
    //                        /|\
    //                       / | \
    //                      /  |  \
    //                     /   |   \
    //                    c  i | j  d
    //                     \   |   /
    //                      \  |  /
    //                       \ | /
    //                        \|/
    //                         a

    OpenMesh::SmartEdgeHandle shared_edge;

    for (auto eh_i : fh_i.edges()) {
        for (auto eh_j : fh_j.edges()) {
            if (eh_i.idx() == eh_j.idx()) {
                shared_edge = eh_i;
            }
        }
    } 

    OpenMesh::SmartVertexHandle vh_a;
    OpenMesh::SmartVertexHandle vh_b;
    OpenMesh::SmartVertexHandle vh_c;
    OpenMesh::SmartVertexHandle vh_d;

    if (shared_edge.v0().idx() < shared_edge.v1().idx()) {
        vh_a = shared_edge.v0();
        vh_b = shared_edge.v1();
    } else {
        vh_a = shared_edge.v1();
        vh_b = shared_edge.v0();
    }
    
    if (shared_edge.halfedge(0).face().idx() == fh_i.idx()) {
        vh_c = shared_edge.halfedge(0).next().to();
        vh_d = shared_edge.halfedge(1).next().to();
    } else {
        vh_c = shared_edge.halfedge(1).next().to();
        vh_d = shared_edge.halfedge(0).next().to();
    }

    ACG::Matrix3x3d Ei = orthogonalize(
        mesh_->point(vh_b) - mesh_->point(vh_a),
        mesh_->point(vh_c) - mesh_->point(vh_a)
    );

    ACG::Matrix3x3d Ej = orthogonalize(
        mesh_->point(vh_b) - mesh_->point(vh_a),
        mesh_->point(vh_b) - mesh_->point(vh_d)
    );

    double angle = mesh_->property(edge_adjuestment_angle_, shared_edge);

    ACG::Matrix3x3d R({
        cos(angle), -sin(angle), 0.,
        sin(angle), cos(angle), 0.,
        0., 0., 1.
    });

    return Ej * R * Ei.inverse() * w0;
}


void VectorFields::initializePlugin() {
    QWidget* toolBox = new QWidget();

    QPushButton* initMeshButton = new QPushButton("Initialize Mesh", toolBox);

    QLabel* singularityLabel = new QLabel("Singularities", toolBox);

    singularityTable = new QTableWidget(0, 2, toolBox);
    singularityTable->setHorizontalHeaderItem(0, new QTableWidgetItem("Vertex id"));
    singularityTable->setHorizontalHeaderItem(1, new QTableWidgetItem("Index"));
    singularityTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    singularityTable->setSelectionMode(QAbstractItemView::SingleSelection);

    QPushButton* addSingularityButton = new QPushButton("Add", toolBox);
    QPushButton* removeSingularityButton = new QPushButton("Remove", toolBox);
    QPushButton* clearSingularityButton = new QPushButton("Clear", toolBox);

    QLabel* constraintsLabel = new QLabel("Directional Constraints", toolBox);
    QPushButton* setConstrainedFacesButton = new QPushButton("Set Constrained Faces", toolBox);
    QLabel* faceVectorLabel = new QLabel("Face Vector", toolBox);
    QLineEdit* setFaceVectorLineEdit = new QLineEdit("0.0 1.0 0.0", toolBox);
    setFaceVectorLineEdit->setPlaceholderText("v0 v1 v2");

    QLabel* trivialConnectionLabel = new QLabel("Trivial Connections", toolBox);
    QPushButton* runButton = new QPushButton("Run!", toolBox);

    QPushButton* showVectorFieldButton = new QPushButton("Show Vector Field", toolBox);

    QPushButton* fakeVectorFieldButton = new QPushButton("Generate Fake Vector Field", toolBox);

    QGridLayout* layout = new QGridLayout(toolBox);
    layout->addWidget(initMeshButton, 0, 1, 1, 2);
    layout->addWidget(singularityLabel, 0, 0, 1, 1);
    layout->addWidget(singularityTable, 1, 0, 1, 3);

    layout->addWidget(addSingularityButton, 2, 0);
    layout->addWidget(removeSingularityButton, 2, 1);
    layout->addWidget(clearSingularityButton, 2, 2);

    layout->addWidget(constraintsLabel, 3, 0, 1, 3);
    layout->addWidget(setConstrainedFacesButton, 4, 0, 1, 3);
    layout->addWidget(faceVectorLabel, 5, 0, 1, 1);
    layout->addWidget(setFaceVectorLineEdit, 5, 1, 1, 2);

    layout->addWidget(trivialConnectionLabel, 6, 0, 1, 3);
    layout->addWidget(runButton, 7, 0, 1, 3);
    layout->addWidget(showVectorFieldButton, 8, 0, 1, 3);

    layout->addWidget(fakeVectorFieldButton, 9, 0, 1, 3);


    // connect(smoothButton, SIGNAL(clicked()), this, SLOT(simpleLaplace()));

    connect(initMeshButton, SIGNAL(clicked()), this, SLOT(initializeMesh()));
    connect(addSingularityButton, SIGNAL(clicked()), this, SLOT(addSingularity()));
    connect(removeSingularityButton, SIGNAL(clicked()), this, SLOT(removeSingularity()));
    connect(clearSingularityButton, SIGNAL(clicked()), this, SLOT(clearSingularities()));

    connect(singularityTable->selectionModel(), SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)), this, SLOT(onItemSelected(const QItemSelection&, const QItemSelection&)));

    connect(singularityTable, SIGNAL(cellChanged(int, int)), this, SLOT(onCellChanged(int, int)));


    connect(runButton, SIGNAL(clicked()), this, SLOT(runAll()));

    connect(showVectorFieldButton, SIGNAL(clicked()), this, SLOT(showVectorField()));
    connect(fakeVectorFieldButton, SIGNAL(clicked()), this, SLOT(generateFakeVectorField()));

    emit addToolbox(tr("VectorField"), toolBox);

}

void VectorFields::initializeMesh() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto trimesh = tri_obj->mesh();

        if (trimesh) {

            mesh_ = trimesh;

            tri_obj->materialNode()->set_point_size(12);

            if(!mesh_->get_property_handle(edge_adjuestment_angle_, "edge adjustment angle")) {
                mesh_->add_property(edge_adjuestment_angle_, "edge adjustment angle");
            }
            for(auto eh : mesh_->edges())
                mesh_->property(edge_adjuestment_angle_, eh) = 0.;

            if(!mesh_->get_property_handle(face_tangent_vector_, "face tangent vector")) {
                mesh_->add_property(face_tangent_vector_, "face tangent vector");
            }
            for(auto fh : mesh_->faces())
                mesh_->property(face_tangent_vector_, fh) = ACG::Vec3d(0, 0, 0);

            
            if(!mesh_->get_property_handle(parent_prime_, "parent prime")) {
                mesh_->add_property(parent_prime_, "parent prime");
            }

            if(!mesh_->get_property_handle(parent_dual_, "parent dual")) {
                mesh_->add_property(parent_dual_, "parent dual");
            }


            tri_obj->meshNode()->drawMode(ACG::SceneGraph::DrawModes::WIREFRAME | ACG::SceneGraph::DrawModes::SOLID_SMOOTH_SHADED | ACG::SceneGraph::DrawModes::POINTS_COLORED);

            tri_obj->materialNode()->enable_alpha_test(0.8);

            updateVertexColors(*mesh_);

            emit updatedObject(tri_obj->id(), UPDATE_ALL);

            std::cout << "Initialized!" << std::endl;
        }
    }


}

void VectorFields::updateVertexColors(TriMesh & _mesh) {


    //reset
    for(auto vh : _mesh.vertices()) {
        if(_mesh.color(vh) != ACG::Vec4f(1,1,0,1)) {
            _mesh.set_color(vh, ACG::Vec4f(1,1,1,0));
        }
    }

    for (int i=0; i<singularities_.size(); i++){
        _mesh.set_color(_mesh.vertex_handle(singularities_[i].first), ACG::Vec4f(0,0,1,1));
    }
}


void VectorFields::refreshSingularityTable() {
    enableOnCellChanged = false;

    singularityTable->setRowCount(0);

    for (int i=0; i<singularities_.size(); i++) {
        singularityTable->insertRow(i);

        QTableWidgetItem* vertex_id_item = new QTableWidgetItem(std::to_string(singularities_[i].first).c_str());
        vertex_id_item->setFlags(vertex_id_item->flags() & (~Qt::ItemIsEditable));
        singularityTable->setItem(i, 0, vertex_id_item);
        singularityTable->setItem(i, 1, new QTableWidgetItem(std::to_string(singularities_[i].second).c_str()));
    }

    enableOnCellChanged = true;

}

void VectorFields::addSingularity() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto trimesh = tri_obj->mesh();

        if (trimesh) {

            // mesh_ = trimesh;

            tri_obj->materialNode()->set_point_size(12);

            for(auto vh : trimesh->vertices()) {
                if(trimesh->status(vh).selected()) {
                    
                    bool skip = false;
                    for (int i=0; i<singularities_.size(); i++) {
                        if (singularities_[i].first == vh.idx()) {
                            skip = true;
                            break;
                        }
                    }
                    if (skip) {
                        continue;
                    }
                    singularities_.push_back(std::pair<int, float> { vh.idx(), 0 });
                    trimesh->status(vh).set_selected(false);
                }
            }

            tri_obj->meshNode()->drawMode(ACG::SceneGraph::DrawModes::WIREFRAME | ACG::SceneGraph::DrawModes::SOLID_SMOOTH_SHADED | ACG::SceneGraph::DrawModes::POINTS_COLORED);

            tri_obj->materialNode()->enable_alpha_test(0.8);

            updateVertexColors(*trimesh);

            emit updatedObject(tri_obj->id(), UPDATE_ALL);

        }
    }

    refreshSingularityTable();

}

void VectorFields::removeSingularity() {

    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto trimesh = tri_obj->mesh();

        if (trimesh) {

            mesh_ = trimesh;

            tri_obj->materialNode()->set_point_size(12);

            for(auto vh : trimesh->vertices()) {
                if(trimesh->status(vh).selected()) {
                    
                    bool skip = true;
                    int index = -1;
                    for (int i=0; i<singularities_.size(); i++) {
                        if (singularities_[i].first == vh.idx()) {
                            skip = false;
                            index = i;
                            break;
                        }
                    }
                    if (!skip) {
                        singularities_.erase(singularities_.begin() + index);
                    }
                    trimesh->status(vh).set_selected(false);
                }
            }

            updateVertexColors(*trimesh);
            emit updatedObject(tri_obj->id(), UPDATE_ALL);

        }
    }

    refreshSingularityTable();

}

void VectorFields::clearSingularities() {
    singularities_.clear();

    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto trimesh = tri_obj->mesh();
        if (trimesh) {
            updateVertexColors(*trimesh);
        }

        emit updatedObject(tri_obj->id(), UPDATE_ALL);
    }

    refreshSingularityTable();

}

void VectorFields::onItemSelected(const QItemSelection& selected, const QItemSelection& unselected) {

    if (selected.size() > 0) {
        // Only 1 item allowed per selection, so we just pick the first one.

        int index = selected[0].top();

        for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {

            auto tri_obj = PluginFunctions::triMeshObject(*o_it);
            auto trimesh = tri_obj->mesh();
            
            if (trimesh) {
                //reset
                for(auto vh : mesh_->vertices()) {
                    mesh_->status(vh).set_selected(false);
                }
                mesh_->status(mesh_->vertex_handle(singularities_[index].first)).set_selected(true);

            }

            emit updatedObject(tri_obj->id(), UPDATE_ALL);
        }

    }
}

void VectorFields::onCellChanged(int row, int column) {

    if (!enableOnCellChanged) {
        return;
    }

    QString text = singularityTable -> item(row, column) -> text();

    bool isFloat = false;
    float new_index = text.toFloat(&isFloat);
    if (isFloat) {
        singularities_[row].second = new_index;
    }

    refreshSingularityTable();

}


////// Jordan's Part

bool VectorFields::inPrimalSpanningTree(TriMesh& mesh_, OpenMesh::HalfedgeHandle he) {
    OpenMesh::VertexHandle v = mesh_.from_vertex_handle(he);
    OpenMesh::VertexHandle w = mesh_.from_vertex_handle(mesh_.opposite_halfedge_handle(he));
    return mesh_.property(parent_prime_, v) == w || mesh_.property(parent_prime_, w) == v;
}

bool VectorFields::inDualSpanningTree(TriMesh& mesh_, OpenMesh::HalfedgeHandle he) {
    OpenMesh::FaceHandle f = mesh_.face_handle(he);
    OpenMesh::FaceHandle g = mesh_.face_handle(mesh_.opposite_halfedge_handle(he));
    return mesh_.property(parent_dual_, f) == g || mesh_.property(parent_dual_, g) == f;
}

void VectorFields::buildPrimalSpanningTree(TriMesh& mesh_) {
    //set root to non-boundary vertex
    auto v_it_root = mesh_.vertices_begin();
    while ((*v_it_root).is_boundary()) v_it_root++;
    OpenMesh::VertexHandle root = *v_it_root;

    //initialize each vertex's parent to itself
        for(auto vh : mesh_.vertices()) {
            mesh_.property(parent_prime_, vh) = vh;
        }
        std::queue<OpenMesh::VertexHandle> Q;
        Q.push(*v_it_root);
        while(!Q.empty()) {
            OpenMesh::VertexHandle v = Q.front(); Q.pop();
            OpenMesh::HalfedgeHandle he =*(mesh_.voh_iter(v)); //get an outgoing halfedge
            OpenMesh::HalfedgeHandle v_ohe = he;
            do {
                //get the originating vertex from the opposite halfedge
                OpenMesh::VertexHandle w = mesh_.from_vertex_handle(mesh_.opposite_halfedge_handle(he));

                if(mesh_.property(parent_prime_, w) == w &&
                    w != root &&
                    !(mesh_.is_boundary(w))) {
                        mesh_.property(parent_prime_, w) = v;
                        Q.push(w);
                }
                //get the next halfedge of the opposite halfedge
                he = mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(he));
            } while (he != v_ohe);
        }
}

void VectorFields::buildDualSpanningCoTree(TriMesh& mesh_) {
    auto f_it_root = mesh_.faces_begin();
    while ((*f_it_root).is_boundary()) f_it_root++;
    OpenMesh::FaceHandle root = *f_it_root;

    for(auto fh : mesh_.faces()) {
        mesh_.property(parent_dual_, fh) = fh;
    }
    std::queue<OpenMesh::FaceHandle> Q;
    Q.push(root);
    while(!Q.empty()) {
        OpenMesh::FaceHandle f = Q.front(); Q.pop();
        OpenMesh::HalfedgeHandle he = *(mesh_.fh_iter(f));
        OpenMesh::HalfedgeHandle f_ohe = he;
        do {
            OpenMesh::FaceHandle g = mesh_.face_handle(mesh_.opposite_halfedge_handle(he));
            if (mesh_.property(parent_dual_, g) == g &&
            g != root &&
            !mesh_.is_boundary(g) && !inPrimalSpanningTree(mesh_, he)) {
                mesh_.property(parent_dual_, g) = f;
                Q.push(g);
            }
        he = mesh_.next_halfedge_handle(he);
        } while (he != f_ohe);
    }
}

void VectorFields::buildTreeCotreeDecomposition(TriMesh& mesh_) {
    buildPrimalSpanningTree(mesh_);
    buildDualSpanningCoTree(mesh_);
}

OpenMesh::HalfedgeHandle VectorFields::sharedHalfEdge(TriMesh& mesh_, OpenMesh::VertexHandle v, OpenMesh::VertexHandle w) {
    OpenMesh::HalfedgeHandle heh = *(mesh_.voh_iter(v));
    OpenMesh::HalfedgeHandle starting_heh = heh;
    do {
        if(mesh_.from_vertex_handle(mesh_.opposite_halfedge_handle(heh)) == w) {
            return heh;
        }
        heh = mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(heh));
    } while (heh != starting_heh);
    assert( 0 );
}

OpenMesh::HalfedgeHandle VectorFields::sharedHalfEdge(TriMesh& mesh_, OpenMesh::FaceHandle f, OpenMesh::FaceHandle g) {
    OpenMesh::HalfedgeHandle heh = *mesh_.fh_begin(f);
    OpenMesh::HalfedgeHandle starting_heh = heh;
    do {
        if(mesh_.face_handle(mesh_.opposite_halfedge_handle(heh)).idx() == g.idx()) {
            return heh;
        }
        heh = mesh_.next_halfedge_handle(heh);
    } while (heh != starting_heh);
    assert( 0 );
}

bool VectorFields::isDualBoundaryLoop(TriMesh& mesh_, const Cycle& cycle) {
    if(cycle.size() == 0) return false;
    return mesh_.is_boundary(mesh_.from_vertex_handle(cycle[0])) || 
        mesh_.is_boundary(mesh_.from_vertex_handle(mesh_.opposite_halfedge_handle(cycle[0])));
}

void VectorFields::appendDualGenerators(TriMesh& mesh_, std::vector<Cycle>& cycles) {
    int orphan_count = 0;
    for(auto eh : mesh_.edges()) {
        if(mesh_.is_boundary(eh)) continue;
        OpenMesh::HalfedgeHandle heh = mesh_.halfedge_handle(eh, 0);
        if(!inPrimalSpanningTree(mesh_, heh) && !inDualSpanningTree(mesh_, heh)) {
            orphan_count++;
            std::cout << "Find " << orphan_count << " edges." << std::endl;
            Cycle g, c1, c2;
            OpenMesh::FaceHandle f;
            g.push_back(heh);
            f = mesh_.face_handle(mesh_.opposite_halfedge_handle(heh));
            while (f != mesh_.property(parent_dual_, f)) {
                c2.push_back(sharedHalfEdge(mesh_, f, mesh_.property(parent_dual_, f)));
                f = mesh_.property(parent_dual_, f);
            }
            f = mesh_.face_handle(heh);
            while (f != mesh_.property(parent_dual_, f)) {
                c1.push_back(sharedHalfEdge(mesh_, f, mesh_.property(parent_dual_, f)));
                f = mesh_.property(parent_dual_, f);
            }
            int m = c1.size()-1;
            int n = c2.size()-1;
            while(c1[m] == c2[n]) {
                m--;
                n--;
            }
            for(int i = 0; i <= m; i++) {
                g.push_back( c1[i] );
            }
            for(int i = n; i >= 0; i--) {
                g.push_back( mesh_.opposite_halfedge_handle(c2[i]));
            }
            //ensure that boundary loops wind around the boundary in a consistent direction
            if (isDualBoundaryLoop(mesh_, g)) {
                if(mesh_.is_boundary(mesh_.from_vertex_handle(mesh_.next_halfedge_handle(g[0])))) {
                    unsigned int n = g.size();
                    for(unsigned int i = 0; i < n; i++) {
                        g[i] = mesh_.opposite_halfedge_handle(g[i]);
                    }
                    for(unsigned int i = 0; i < n/2; i++) {
                        std::swap(g[i], g[n-1-i]);
                    }
                }
            }
            cycles.push_back(g);
        }

    }
}

////// Jordan's Part...END


void VectorFields::runAll() {
    // Pipeline:
    // 1. Find contractible & non- cycles
    // 2. Compute angle defects for each cycle
    // 3. Build Matrix A=[d0 H]^T and vector b
    // 4. Replace rows in A and b for singularities
    // (Extra) Add rows for directional constraints
    // 5. Build the linear system, and solve for adjustment angles
    // 6. Construct the direction fields, save the angle with each face
    // 7. Visualise the vector field on the mesh

    std::vector<Cycle> non_contra_cycles;
    buildTreeCotreeDecomposition(*mesh_);

    // showPrimalTree();
    // showDualTree();

    appendDualGenerators(*mesh_, non_contra_cycles);


    std::cout << non_contra_cycles.size() << std::endl;
    // std::cout << non_contra_cycles[0].size() << std::endl;
    // std::cout << non_contra_cycles[10].size() << std::endl;

    showCycles(non_contra_cycles);

}

void VectorFields::showPrimalTree() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        if (tri_obj->mesh() == mesh_) {

            ACG::SceneGraph::LineNode* lineNode;

            if (!tri_obj->getAdditionalNode(lineNode, name(),"Primal Tree"))
            {
                lineNode = new ACG::SceneGraph::LineNode(ACG::SceneGraph::LineNode::LineSegmentsMode,
                        tri_obj->manipulatorNode(),"Primal Tree");
                tri_obj->addAdditionalNode(lineNode, name(), "Primal Tree");

                //creates the line
                lineNode->clear_points();
                lineNode->set_color(OpenMesh::Vec4f(1.0f,0.0f,0.0f,1.0f));
                lineNode->set_line_width(3);
                // lineNode->add_line(p0, p1);
                // lineNode->alwaysOnTop() = true;
            }else {
                //creates the line
                lineNode->clear_points();
                lineNode->set_color(OpenMesh::Vec4f(1.0f,0.0f,0.0f,1.0f));
                lineNode->set_line_width(3);
                // lineNode->add_line(p0, p1);
                // lineNode->alwaysOnTop() = true;
            }

            for (auto vh : mesh_->vertices()) {
                if (mesh_->property(parent_prime_, vh) != vh) {

                    lineNode->add_color(OpenMesh::Vec4f(0.0f,0.0f,1.0f,1.0f));
                    lineNode->add_line(mesh_->point(vh), mesh_->point(mesh_->property(parent_prime_, vh)));
                }
            }
            
        }
    }

    emit updateView();

}


void VectorFields::showDualTree() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        if (tri_obj->mesh() == mesh_) {

            ACG::SceneGraph::LineNode* lineNode;

            if (!tri_obj->getAdditionalNode(lineNode, name(),"Dual Tree"))
            {
                lineNode = new ACG::SceneGraph::LineNode(ACG::SceneGraph::LineNode::LineSegmentsMode,
                        tri_obj->manipulatorNode(),"Dual Tree");
                tri_obj->addAdditionalNode(lineNode, name(), "Dual Tree");

                //creates the line
                lineNode->clear_points();
                lineNode->set_color(OpenMesh::Vec4f(1.0f,0.0f,0.0f,1.0f));
                lineNode->set_line_width(3);
                // lineNode->add_line(p0, p1);
                // lineNode->alwaysOnTop() = true;
            }else {
                //creates the line
                lineNode->clear_points();
                lineNode->set_color(OpenMesh::Vec4f(1.0f,0.0f,0.0f,1.0f));
                lineNode->set_line_width(3);
                // lineNode->add_line(p0, p1);
                // lineNode->alwaysOnTop() = true;
            }

            for (auto fh : mesh_->faces()) {
                if (mesh_->property(parent_dual_, fh) != fh) {
                    ACG::Vec3d centroid1(0, 0, 0);
                    for (auto vh : fh.vertices()) {
                        centroid1 += mesh_->point(vh);
                    }
                    centroid1 /= 3;

                    ACG::Vec3d centroid2(0, 0, 0);
                    for (auto vh : OpenMesh::SmartFaceHandle(mesh_->property(parent_dual_, fh).idx(), mesh_).vertices()) {
                        centroid2 += mesh_->point(vh);
                    }
                    centroid2 /= 3;

                    lineNode->add_color(OpenMesh::Vec4f(1.0f,0.0f,0.0f,1.0f));
                    lineNode->add_line(centroid1, centroid2);

                    lineNode->add_color(OpenMesh::Vec4f(1.0f,1.0f,0.0f,1.0f));
                    lineNode->add_line(0.8 * centroid2 + 0.2 * centroid1, centroid2);

                }
            }
            
        }
    }

    emit updateView();

}


void VectorFields::showCycles(std::vector<Cycle> & cycles) {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        if (tri_obj->mesh() == mesh_) {

            ACG::SceneGraph::LineNode* lineNode;

            if (!tri_obj->getAdditionalNode(lineNode, name(),"Cycles"))
            {
                lineNode = new ACG::SceneGraph::LineNode(ACG::SceneGraph::LineNode::LineSegmentsMode,
                        tri_obj->manipulatorNode(),"Cycles");
                tri_obj->addAdditionalNode(lineNode, name(), "Cycles");

                //creates the line
                lineNode->clear_points();
                lineNode->set_color(OpenMesh::Vec4f(1.0f,0.0f,0.0f,1.0f));
                lineNode->set_line_width(3);
                // lineNode->add_line(p0, p1);
                // lineNode->alwaysOnTop() = true;
            }else {
                //creates the line
                lineNode->clear_points();
                lineNode->set_color(OpenMesh::Vec4f(1.0f,0.0f,0.0f,1.0f));
                lineNode->set_line_width(3);
                // lineNode->add_line(p0, p1);
                // lineNode->alwaysOnTop() = true;
            }

            for (int i=0; i<cycles.size(); i++) {
            // for (int i=0; i<1; i++) {
                OpenMesh::Vec4f color((double) rand() / RAND_MAX, (double) rand() / RAND_MAX, (double) rand() / RAND_MAX, 1.0);
                for (int j=0; j<cycles[i].size(); j++) {
                // for (int j=0; j<5; j++) {
                    OpenMesh::SmartHalfedgeHandle heh = OpenMesh::SmartHalfedgeHandle(cycles[i][j].idx(), mesh_);

                    OpenMesh::SmartFaceHandle fh1 = heh.face();
                    OpenMesh::SmartFaceHandle fh2 = heh.opp().face();

                    ACG::Vec3d centroid1(0, 0, 0);
                    for (auto vh : fh1.vertices()) {
                        centroid1 += mesh_->point(vh);
                    }
                    centroid1 /= 3;

                    ACG::Vec3d centroid2(0, 0, 0);
                    for (auto vh : fh2.vertices()) {
                        centroid2 += mesh_->point(vh);
                    }
                    centroid2 /= 3;

                    if (j==0) {
                        lineNode->add_color(OpenMesh::Vec4f(0.0f,1.0f,1.0f,1.0f));
                        lineNode->add_line(centroid1, centroid2);
                    } else {
                        lineNode->add_color(color);
                        lineNode->add_line(centroid1, centroid2);
                    }

                    lineNode->add_color(OpenMesh::Vec4f(1.0f,0.0f,1.0f,1.0f));
                    lineNode->add_line(0.8 * centroid2 + 0.2 * centroid1, centroid2);

                }
            }
        }
    }

    emit updateView();
}

void VectorFields::showVectorField() {

    if (!mesh_) return;

    OpenMesh::SmartFaceHandle root_face = *(mesh_->faces_begin());
    OpenMesh::SmartHalfedgeHandle root_heh;

    for (auto heh : root_face.halfedges()) {
        root_heh = heh;
        break;
    }
    ACG::Vec3d init_direction = (mesh_->point(root_heh.to()) - mesh_->point(root_heh.from())).normalize();
    mesh_->property(face_tangent_vector_, root_face) = init_direction;

    std::vector<OpenMesh::SmartFaceHandle> Q;
    std::map<int, bool> visited;
    Q.push_back(root_face);
    visited[root_face.idx()] = true;

    while (Q.size() > 0) {
        OpenMesh::SmartFaceHandle fh_i = Q[Q.size()-1];
        Q.pop_back();

        for (auto fh_j : fh_i.faces()) {
            if (!visited[fh_j.idx()]) {
                ACG::Vec3d dir_j = transport(mesh_->property(face_tangent_vector_, fh_i), fh_i, fh_j);
                mesh_->property(face_tangent_vector_, fh_j) = dir_j;

                Q.insert(Q.begin(), fh_j);
                visited[fh_j.idx()] = true;
            }
        }
    }

    std::cout << "finished computation!" << std::endl;

    // Draw lines over faces:
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        if (tri_obj->mesh() == mesh_) {


            ACG::SceneGraph::LineNode* lineNode;
            //create line node
            if (!tri_obj->getAdditionalNode(lineNode, name(),"Vector Field"))
            {
                lineNode = new ACG::SceneGraph::LineNode(ACG::SceneGraph::LineNode::LineSegmentsMode,
                        tri_obj->manipulatorNode(),"Vector Field");
                tri_obj->addAdditionalNode(lineNode, name(), "Vector Field");

                //creates the line
                lineNode->clear_points();
                lineNode->set_color(OpenMesh::Vec4f(1.0f,0.0f,0.0f,1.0f));
                lineNode->set_line_width(3);
                // lineNode->add_line(p0, p1);
                // lineNode->alwaysOnTop() = true;
            }else {
                //creates the line
                lineNode->clear_points();
                lineNode->set_color(OpenMesh::Vec4f(1.0f,0.0f,0.0f,1.0f));
                lineNode->set_line_width(3);
                // lineNode->add_line(p0, p1);
                // lineNode->alwaysOnTop() = true;
            }

            for (auto fh : mesh_->faces()) {

                double scale = 1;
                for (auto heh : root_face.halfedges()) {
                    scale = (mesh_->point(heh.to()) - mesh_->point(heh.from())).norm() / 2;
                    break;
                }


                ACG::Vec3d dir = mesh_->property(face_tangent_vector_, fh);

                ACG::Vec3d centroid(0, 0, 0);
                for (auto vh : fh.vertices()) {
                    centroid += mesh_->point(vh);
                }
                centroid /= 3;

                ACG::Vec3d p0 = centroid - scale / 2 * dir;
                ACG::Vec3d p1 = centroid + scale / 2 * dir;

                lineNode->add_line(p0, p1);
            }
        }
    }

    emit updateView();
}

void VectorFields::generateFakeVectorField() {
    if (!mesh_) return;

    double increment = 0.;
    double angle = 0;
    for (auto eh : mesh_->edges()) {
        mesh_->property(edge_adjuestment_angle_, eh) = angle * M_PI;
        angle = std::fmod((angle + increment), 2.);
    }

    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {

        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        emit updatedObject(tri_obj->id(), UPDATE_ALL);
    }
    
    std::cout << "finished fake!" << std::endl;

}

