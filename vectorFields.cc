#include "vectorFields.hh"
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"

#include <string>
#include <queue>
#include <algorithm>
#include <chrono>

void showMessage(std::string message) {
    QMessageBox* msgBox = new QMessageBox();

    msgBox->setText(QString::fromStdString(message));
    msgBox->exec();
    return;
}

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

    // if (shared_edge.v0().idx() < shared_edge.v1().idx()) {
    //     vh_a = shared_edge.v0();
    //     vh_b = shared_edge.v1();
    // } else {
    //     vh_a = shared_edge.v1();
    //     vh_b = shared_edge.v0();
    // }
    

    if (shared_edge.h0().face().idx() == fh_i.idx()) {

        vh_a = shared_edge.h0().from();
        vh_b = shared_edge.h0().to();

        vh_c = shared_edge.halfedge(0).next().to();
        vh_d = shared_edge.halfedge(1).next().to();

    } else {

        vh_a = shared_edge.h0().to();
        vh_b = shared_edge.h0().from();

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

    // if (fh_i.idx() > fh_j.idx()) {
    //     angle = - angle;
    // }

    if (shared_edge.h0().face().idx() == fh_i.idx()) {
        angle = - angle;
    }

    // if (shared_edge.idx() == 1080 || shared_edge.idx() == 1081) {
    //     std::cout << (shared_edge.h0().face().idx() == fh_j.idx()) << std::endl;
    //     std::cout << (d1_.coeffRef(fh_i.idx(), shared_edge.idx())) << std::endl;
    //     std::cout << "Edge " << shared_edge.idx() << " angle " << angle << std::endl;
    //     std::cout << "Ei " << Ei << std::endl;
    //     std::cout << "Ej " << Ej << std::endl;
    // }

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
    QLabel* rootFaceAngleLabel = new QLabel("Root face angle", toolBox);

    setRootFaceAngleLineEdit->setPlaceholderText("degrees");

    QLabel* trivialConnectionLabel = new QLabel("Trivial Connections", toolBox);
    QPushButton* runButton = new QPushButton("Run!", toolBox);

    QPushButton* showVectorFieldButton = new QPushButton("Show Vector Field", toolBox);

    QPushButton* fakeVectorFieldButton = new QPushButton("Generate Fake Vector Field", toolBox);

    QGridLayout* layout = new QGridLayout(toolBox);
    layout->addWidget(initMeshButton, 0, 0, 1, 3);
    layout->addWidget(singularityLabel, 1, 0, 1, 3);
    layout->addWidget(singularityTable, 2, 0, 1, 3);

    layout->addWidget(addSingularityButton, 3, 0);
    layout->addWidget(removeSingularityButton, 3, 1);
    layout->addWidget(clearSingularityButton, 3, 2);

    layout->addWidget(constraintsLabel, 4, 0, 1, 3);
    layout->addWidget(setConstrainedFacesButton, 5, 0, 1, 3);
    layout->addWidget(rootFaceAngleLabel, 6, 0, 1, 1);
    layout->addWidget(setRootFaceAngleLineEdit, 6, 1, 1, 2);

    layout->addWidget(trivialConnectionLabel, 7, 0, 1, 3);
    layout->addWidget(runButton, 8, 0, 1, 3);
    layout->addWidget(showVectorFieldButton, 9, 0, 1, 3);

    layout->addWidget(fakeVectorFieldButton, 10, 0, 1, 3);


    // connect(smoothButton, SIGNAL(clicked()), this, SLOT(simpleLaplace()));

    connect(initMeshButton, SIGNAL(clicked()), this, SLOT(initializeMesh()));
    connect(addSingularityButton, SIGNAL(clicked()), this, SLOT(addSingularity()));
    connect(removeSingularityButton, SIGNAL(clicked()), this, SLOT(removeSingularity()));
    connect(clearSingularityButton, SIGNAL(clicked()), this, SLOT(clearSingularities()));

    connect(singularityTable->selectionModel(), SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)), this, SLOT(onItemSelected(const QItemSelection&, const QItemSelection&)));

    connect(singularityTable, SIGNAL(cellChanged(int, int)), this, SLOT(onCellChanged(int, int)));

    connect(setConstrainedFacesButton, SIGNAL(clicked()), this, SLOT(setConstraintFaces()));
    connect(setRootFaceAngleLineEdit, SIGNAL(textEdited(const QString &)), this, SLOT(onRootFaceAngleChanged(const QString &)));

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


            tri_obj->meshNode()->drawMode(ACG::SceneGraph::DrawModes::WIREFRAME | ACG::SceneGraph::DrawModes::SOLID_SMOOTH_SHADED | ACG::SceneGraph::DrawModes::POINTS_COLORED | ACG::SceneGraph::DrawModes::SOLID_FACES_COLORED);

            tri_obj->materialNode()->enable_alpha_test(0.8);
            updateVertexColors(*mesh_);

            emit updatedObject(tri_obj->id(), UPDATE_ALL);

            findCyclesAndBuildA(nullptr);

            constrained_faces_.clear();
            // showMessage("Init complete! \n" );

        }

    }


}

void VectorFields::findCyclesAndBuildA(OpenMesh::SmartFaceHandle* root_fh) {

    uint64_t start_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();


    // Find Cycles

    std::vector<Cycle> contraCycles;
    std::vector<Cycle> nonContraCycles;

    std::cout << "Finding contractible cycles..." << std::endl;
    findContractibleLoops(*mesh_, contraCycles);

    std::cout << "Finding non-contractible cycles..." << std::endl;
    buildTreeCotreeDecomposition(*mesh_, root_fh);
    appendDualGenerators(*mesh_, nonContraCycles);

    nGenerators = nonContraCycles.size();

    std::cout << "Found " << nGenerators << " non-contractible cycles..." << std::endl;


    std::vector<Cycle> constrainedFacePaths;
    findConstrainedFacePaths(*mesh_, constrainedFacePaths);

    std::cout << "Num of paths: " << constrainedFacePaths.size() << std::endl;

    basisCycles.clear();
    basisCycles.reserve(contraCycles.size() + nonContraCycles.size());
    basisCycles.insert(basisCycles.end(), contraCycles.begin(), contraCycles.end());
    basisCycles.insert(basisCycles.end(), nonContraCycles.begin(), nonContraCycles.end());
    basisCycles.insert(basisCycles.end(), constrainedFacePaths.begin(), constrainedFacePaths.end());

    showCycles(nonContraCycles);

    // Build Matrix A
    std::cout << "Building Matrix A..." << std::endl;

    A = Eigen::SparseMatrix<double>(static_cast<size_t>(mesh_->n_edges()), basisCycles.size());
    buildCycleMatrix(A, basisCycles);

    // Compute the defects

    K = Eigen::MatrixXd(Eigen::MatrixXd::Zero(basisCycles.size(), 1));
    b = Eigen::MatrixXd(Eigen::MatrixXd::Zero(basisCycles.size(), 1));
    std::vector<Eigen::Triplet<double>> triplets_K;
    std::vector<Eigen::Triplet<double>> triplets_b;
    // generatorOnBoundary.resize(nGenerators());

    double vert_defect;
    int j = 0;
    for (auto vh : mesh_->vertices()) {
        if (j < contraCycles.size()) {
            vert_defect = vertex_defect(vh);
            K(j) = vert_defect;
        }
        else {
            break;
        }
        j++;
    }
    
    for (int i = 0; i < nonContraCycles.size(); i++) {
        // if (isDualBoundaryLoop(mesh_, basisCycles[i])) {
        //     K( i ) = -boundaryLoopCurvature( basisCycles[i] );
        //     generatorOnBoundary[i - nContractibleCycles ] = true;
        // }
        // else {
            K( i + contraCycles.size() ) = -defect(nonContraCycles[i]);
        //     generatorOnBoundary[i-nContractibleCycles] = false;
        // }
    }

    for (int i = 0; i < constrainedFacePaths.size(); i++) {
        // if (isDualBoundaryLoop(mesh_, basisCycles[i])) {
        //     K( i ) = -boundaryLoopCurvature( basisCycles[i] );
        //     generatorOnBoundary[i - nContractibleCycles ] = true;
        // }
        // else {
            double theta = root_face_angle_ / 180 * M_PI;

            for (auto he : constrainedFacePaths[i]) {
                theta = parallelTransport(theta, he);
            }

            K( i + contraCycles.size() + nonContraCycles.size()) = root_face_angle_ / 180 * M_PI - theta;
        //     generatorOnBoundary[i-nContractibleCycles] = false;
        // }
    }



    A = A.transpose();

    std::cout << "QRSolver: Computing A..." << std::endl;
    QRSolver.compute(A);
    std::cout << "QRSolver: Computed A." << std::endl;

    // Build d1
    buildD1(d1_);
    uint64_t end_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    std::cout << "Initialized! Time: " << end_time - start_time << "ms." << std::endl;
    std::cout << "Matrix A: " << A.rows() << " rows " << A.cols() << " cols." << std::endl;
    std::cout << "Matrix K: " << K.rows() << " rows " << K.cols() << " cols." << std::endl;
    std::cout << "Matrix d1: " << d1_.rows() << " rows " << d1_.cols() << " cols." << std::endl;

    showMessage("Found " + std::to_string(contraCycles.size()) + " contra cycles and " + std::to_string(nGenerators) + " non-contra cycles." );

}

void VectorFields::findConstrainedFacePaths(TriMesh & _mesh, std::vector<Cycle> & paths) {
    if (constrained_faces_.size() > 1) {
        OpenMesh::SmartFaceHandle root_face = constrained_faces_[0];
        for (int i=1; i<constrained_faces_.size(); i++) {
            Cycle path;
            OpenMesh::SmartFaceHandle fh = constrained_faces_[i];

            do {
                for (auto heh:fh.halfedges()) {
                    if (heh.opp().face().idx()==mesh_->property(parent_dual_, fh).idx()){
                        path.push_back(heh);
                        fh = heh.opp().face();
                        break;
                    }
                }

            } while (std::find(constrained_faces_.begin(), constrained_faces_.end(), fh) == constrained_faces_.end());


            paths.push_back(path);
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

        if (mesh_ == trimesh) {

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

            updateVertexColors(*trimesh);

            emit updatedObject(tri_obj->id(), UPDATE_ALL);

            refreshSingularityTable();
            return;
        }
    }

    showMessage("Please initialize the mesh first.");

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

void VectorFields::buildDualSpanningCoTree(TriMesh& mesh_, OpenMesh::SmartFaceHandle* root_fh) {
    
    OpenMesh::FaceHandle root;
    if (root_fh) {
        root = *root_fh;
    } else {
        auto f_it_root = mesh_.faces_begin();
        while ((*f_it_root).is_boundary()) f_it_root++;
        root = *f_it_root;
    }



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

void VectorFields::buildTreeCotreeDecomposition(TriMesh& mesh_, OpenMesh::SmartFaceHandle* root_fh) {
    buildPrimalSpanningTree(mesh_);
    buildDualSpanningCoTree(mesh_, root_fh);
}

OpenMesh::SmartHalfedgeHandle VectorFields::sharedHalfEdge(TriMesh& mesh_, OpenMesh::VertexHandle v, OpenMesh::VertexHandle w) {
    OpenMesh::HalfedgeHandle heh = *(mesh_.voh_iter(v));
    OpenMesh::HalfedgeHandle starting_heh = heh;
    do {
        if(mesh_.from_vertex_handle(mesh_.opposite_halfedge_handle(heh)) == w) {
            return OpenMesh::SmartHalfedgeHandle(heh.idx(), &mesh_);
        }
        heh = mesh_.next_halfedge_handle(mesh_.opposite_halfedge_handle(heh));
    } while (heh != starting_heh);
    assert( 0 );
}

OpenMesh::SmartHalfedgeHandle VectorFields::sharedHalfEdge(TriMesh& mesh_, OpenMesh::FaceHandle f, OpenMesh::FaceHandle g) {
    OpenMesh::HalfedgeHandle heh = *mesh_.fh_begin(f);
    OpenMesh::HalfedgeHandle starting_heh = heh;
    do {
        if(mesh_.face_handle(mesh_.opposite_halfedge_handle(heh)).idx() == g.idx()) {
            return OpenMesh::SmartHalfedgeHandle(heh.idx(), &mesh_);
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
            g.push_back(OpenMesh::SmartHalfedgeHandle(heh.idx(), &mesh_));
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


////// Will's Part

// add 2*pi*k to the right hand side, where k is the vector of singularity/generator indices
void VectorFields::setupRHS() {
    double indexSum = 0;

    // iterate over vertices
    for (auto vh : mesh_->vertices()) {

    }
}

void VectorFields::resetRHS() {
    // make a cop
      // make a copy of the right hand side used in the most recent solve
    for( int i = 0; i < basisCycles.size(); i++ )
      {
        b(i) = K(i);
    }
}

void VectorFields::update() {
    //... not sure what to do with Dense x, y?
}

double VectorFields::tipAngle(OpenMesh::Vec3d& x, OpenMesh::Vec3d& a, OpenMesh::Vec3d& b) 
// returns the angle between (a-x) and (b-x)
{
    auto u = (a - x).normalize();
    auto v = (b - x).normalize();

    return atan2( OpenMesh::cross(u,v).norm(), OpenMesh::dot(u, v));
}


double VectorFields::boundaryLoopCurvature(Cycle& cycle) {
    double totalK = 0.;

    // get a halfedge of the "virtual" face bounded by the current cycle
    auto v0 = cycle[0].opp().next().from();
    auto he0 = v0.out();
    do
    {
        he0 = he0.opp().next();
    }
        // is_boundary might not be the same as onBoundary in their implementation...
    while (!he0.is_boundary());

    // compute a "virtual" vertex in the middle of this loop
    OpenMesh::Vec3d c(0., 0., 0.);
    auto he = he0;
    int boundaryLength = 0;
    do
    {
        c += mesh_->point(he.from());
        boundaryLength++;
        he = he.next();
    }
    while (he != he0);
    c /= (double) boundaryLength;

    // compute the curvature around the center vertex
    double K = 2.*M_PI;
    he = he0;
    do
    {
        OpenMesh::Vec3d a = mesh_->point(he.from());
        OpenMesh::Vec3d b = mesh_->point(he.next().from());
        K -= tipAngle(c, a, b);
        he = he.next();
    }
    while (he != he0);
    totalK += K;

    // add the curvature around each of the boundary vertices, using
    // the following labels:
    //    c - virtual center vertex of boundary loop (computed above)
    //    d - current boundary vertex (we walk around the 1-ring of this vertex)
    //    a,b - consecutive interior vertices in 1-ring of d
    //    e,f - boundary vertices adjacent to d
    he = he0;
    do
    {
        auto v = he.from();
        OpenMesh::Vec3d d = mesh_->point(v);

        K = 2.*M_PI;

        auto he2 = v.out();
        do
        {
            if (he2.is_boundary())
            {
                auto f = mesh_->point(he2.next().from());
                K -= tipAngle(d, f, c);
            }
            else
            {
                auto a = mesh_->point(he2.next().from());
                auto b = mesh_->point(he2.next().next().from());
                K -= tipAngle(d, a, b);

                if (he.opp().is_boundary())
                {
                    auto e = mesh_->point(he2.opp().from());
                    K -= tipAngle( d, c, e );
                }
            }

            he2 = he2.opp().next();
        }
        while ( he2 != v.out() );

        totalK += K;

        he = he.next();
    }
    while (he != he0);

    return totalK;
}

double VectorFields::vertex_defect(OpenMesh::SmartVertexHandle vh) {
    double sum = 0.;

    // iterate over incident triangles
    auto he = vh.out();
    // wtf is this while loop structure lmao
    do
    {
        // grab vertices
        auto p1 = mesh_->point(he.from());
        auto p2 = mesh_->point(he.next().from());
        auto p3 = mesh_->point(he.next().next().from());

        // subtract incident angle from sum
        Vector u1 = (p2 - p1);
        Vector u2 = (p3 - p1);
        sum += atan2( (OpenMesh::cross(u1,u2)).norm(), OpenMesh::dot(u1, u2));

        he = he.opp().next();
    }
    while (he != vh.out());
    return 2.*M_PI - sum;
}

double VectorFields::parallelTransport(double phi, OpenMesh::SmartHalfedgeHandle he) 
   // given an angle phi relative to the canonical reference frame
   // of he->face, returns the angle parallel transported across he
   // using the Levi-Civita connection, expressed relative to the
   // canonical frame of he->flip->face
{
    // get (oriented) direction along shared edge;
    auto u = he.from();
    auto v = he.opp().from();
    auto e = mesh_->point(v) - mesh_->point(u);
    if (u.idx() > v.idx()) e = -e;

    // compute angle adjustments between canonical frames;
    ACG::Vec3d a = mesh_->point(he.from());
    ACG::Vec3d b = mesh_->point(he.next().from());
    ACG::Vec3d c = mesh_->point(he.next().next().from());
    OpenMesh::Vec3d e1, e2;
            // is .unit() the same as normalizing?
    e1 = (b - a).normalize();
    e2 = c - a;
    e2 = (e2 - (e2*e1)*e1).normalize();

    he = he.opp();
    a = mesh_->point(he.from());
    b = mesh_->point(he.next().from());
    c = mesh_->point(he.next().next().from());
    OpenMesh::Vec3d f1, f2;
            // is .unit() the same as normalizing?
    f1 = (b - a).normalize();
    f2 = c - a;
    f2 = (f2 - (f2*f1)*f1).normalize();

    double deltaIJ = atan2(OpenMesh::dot(e, e2), OpenMesh::dot(e, e1));
    double deltaJI = atan2(OpenMesh::dot(e, f2), OpenMesh::dot(e, f1));

    // transport phi
    return (phi - deltaIJ) + deltaJI;
}

double VectorFields::defect(Cycle& c) {
    double theta = 0.;

    for (auto he : c) {
        theta = parallelTransport(theta, he);
    }

    while (theta >= M_PI) theta -= 2.*M_PI;
    while (theta < -M_PI) theta += 2.*M_PI;

    return -theta;

}


void VectorFields::buildCycleMatrix(Eigen::SparseMatrix<double>& A, std::vector<Cycle>& cycles) {

    std::vector< Eigen::Triplet<double> > triplets_A;

    for( unsigned int l = 0; l < cycles.size(); l++ )
    {
        for(  auto h  = cycles[l].begin();
                h != cycles[l].end();
                h ++ )
        {
            OpenMesh::SmartHalfedgeHandle heh = *h;
            auto vh = mesh_->to_vertex_handle(heh);
            // these indices might be wrong...
            int k = heh.edge().idx();
            int i = heh.face().idx();
                // might not give the same half edge as in their code...
            int j = heh.opp().face().idx();

            // if( i > j ) {
            //     triplets_A.emplace_back(k, l, -1.);
            // }
            // else {
            //     triplets_A.emplace_back(k, l, 1.);
            // }

            if (heh.idx() == heh.edge().h0().idx()) {
                triplets_A.emplace_back(k, l, 1.);
            }
            else {
                triplets_A.emplace_back(k, l, -1.);
            }
        }
    }
    A.setFromTriplets(triplets_A.begin(), triplets_A.end());
}

void VectorFields::appendDirectionalConstraints(TriMesh& mesh, std::vector<Cycle>& basisCycles, std::vector<double>& holonomies) {
    // first point all faces to themselves to indicate that they have not yet
    // been added to the tree; meanwhile look for a constrained face to serve
    // as the root for our constraint tree (if there aren't any constrained
    // faces, an arbitrary face will work just fine)
    // for ( auto fh : mesh.faces() ) {
        
    // }
}

void VectorFields::findContractibleLoops(TriMesh& mesh, std::vector<Cycle>& basisCycles) {
    // contractible bases
                    // maybe n_vertices = size()?
    vertex2row.resize(mesh.n_vertices());
    for (auto vh : mesh_->vertices()) {
        // std::cout << "Vertex id: " << vh.idx() << std::endl;
        if (vh.is_boundary())
        {
            // maybe subtract one, says theirs is 0-based? not sure what that means
            vertex2row[ vh.idx() ] = -1;
            continue;
        }
        Cycle c;
        // returns an outgoing half-edge
        auto he = vh.halfedge();
        do
        {
            c.push_back( he );
            //   he->flip->next
            he = he.opp().next();
        }
        while ( he != vh.out());

        vertex2row[ vh.idx() ] = static_cast<int>(basisCycles.size()-1);
        basisCycles.push_back(c);
    }
}

////// Will's Part End

void VectorFields::buildD1(Eigen::SparseMatrix<double>& d1) {
    d1.resize(mesh_->n_faces(), mesh_->n_edges());
    std::vector< Eigen::Triplet<double> > triplets_d1;

    for (auto fh : mesh_->faces()) {
        int i = fh.idx();

        auto heh = fh.halfedge();
        int start_heh = heh.idx();
        while (true) {
            int j = heh.opp().face().idx();
            // if (i<j) {
            //     triplets_d1.emplace_back(i, heh.edge().idx(), 1.);
            // } else {
            //     triplets_d1.emplace_back(i, heh.edge().idx(), -1);
            // }

            if (heh.idx() == heh.edge().h0().idx()) {
                triplets_d1.emplace_back(i, heh.edge().idx(), 1.);
            } else {
                triplets_d1.emplace_back(i, heh.edge().idx(), -1);
            }

            heh = heh.next();
            if (heh.idx()==start_heh) {
                break;
            } 
        }

    }

    d1.setFromTriplets(triplets_d1.begin(), triplets_d1.end());
}

void VectorFields::setConstraintFaces() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        auto trimesh = tri_obj->mesh();

        if (mesh_ == trimesh) {

            tri_obj->materialNode()->set_point_size(12);

            constrained_faces_.clear();
            for(auto fh : trimesh->faces()) {
                mesh_->set_color(fh, ACG::Vec4f(1,1,1,0)); // reset

                if(trimesh->status(fh).selected()) {
                    if (std::find(constrained_faces_.begin(), constrained_faces_.end(), fh) != constrained_faces_.end()){
                        continue;
                    }
                    constrained_faces_.push_back(fh);
                    trimesh->status(fh).set_selected(false);
                    mesh_->set_color(fh, ACG::Vec4f(0,1,0,1));

                    mesh_->status(fh).set_selected(false);
                }

            }

            emit updatedObject(tri_obj->id(), UPDATE_ALL);

            std::cout << "Constrained faces: " << constrained_faces_.size() << std::endl;

        }
    }

    if (constrained_faces_.size() > 0) {
        showRootFaceVector();
        findCyclesAndBuildA(&constrained_faces_[0]);
    }

}

void VectorFields::onRootFaceAngleChanged(const QString & text) {
    if (text == "") {
        root_face_angle_ = 0;
        setRootFaceAngleLineEdit->setText("0");
    }

    bool isInt = false;
    float new_angle = text.toInt(&isInt);
    if (isInt) {
        root_face_angle_ = new_angle;
    } else {
        setRootFaceAngleLineEdit->setText(std::to_string(root_face_angle_).c_str());
    }

    showRootFaceVector();
}

void VectorFields::showRootFaceVector() {
    if (constrained_faces_.size() == 0) {
        return;
    }

    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        if (tri_obj->mesh() == mesh_) {


            ACG::SceneGraph::LineNode* lineNode;
            //create line node
            if (!tri_obj->getAdditionalNode(lineNode, name(), "Root Face Vector"))
            {
                lineNode = new ACG::SceneGraph::LineNode(ACG::SceneGraph::LineNode::LineSegmentsMode,
                        tri_obj->manipulatorNode(),"Root Face Vector");
                tri_obj->addAdditionalNode(lineNode, name(), "Root Face Vector");

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

            OpenMesh::SmartFaceHandle root_face = constrained_faces_[0];
            OpenMesh::SmartHalfedgeHandle root_heh = root_face.halfedge();

            double scale = 1;
            for (auto heh : root_face.halfedges()) {
                scale = (mesh_->point(heh.to()) - mesh_->point(heh.from())).norm() / 2;
                break;
            }

            ACG::Vec3d init_direction = (mesh_->point(root_heh.to()) - mesh_->point(root_heh.from())).normalize();

            ACG::Matrix3x3d E = orthogonalize(
                mesh_->point(root_heh.from()) - mesh_->point(root_heh.to()),
                mesh_->point(root_heh.next().to()) - mesh_->point(root_heh.next().from())
            );

            double rad = -root_face_angle_ / 180. * M_PI;
            ACG::Matrix3x3d R({
                cos(rad), -sin(rad), 0.,
                sin(rad), cos(rad), 0.,
                0., 0., 1.
            });



            ACG::Vec3d dir = E * R * E.inverse() * init_direction;

            root_face_dir_ = dir;

            ACG::Vec3d centroid(0, 0, 0);
            for (auto vh : root_face.vertices()) {
                centroid += mesh_->point(vh);
            }
            centroid /= 3;

            ACG::Vec3d p0 = centroid - scale / 2. * dir;
            ACG::Vec3d p1 = centroid + scale / 2. * dir;

            lineNode->add_color(OpenMesh::Vec4f(0.0f,0.0f,1.0f,1.0f));
            lineNode->add_line(p0, p1);

            lineNode->add_color(OpenMesh::Vec4f(0.0f,1.0f,0.0f,1.0f));
            lineNode->add_line(.2*p0+.8*p1, p1);

        }
    }

    emit updateView();

}

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

    uint64_t start_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    bool initialized = false;
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {
        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        if (mesh_ == tri_obj->mesh()) {
            initialized = true;
        }
    }

    if (!initialized) {
        showMessage("Please initialize the mesh first.");
        return;
    }


    // Deal with singularities, build vector b
    // First copy K to b
    std::cout << "Copying K to b..." << std::endl;
    for (int i = 0; i < basisCycles.size(); i++) {
        b(i) = K(i);
    }
    
    // Then update b
    std::cout << "Updating b with singularities..." << std::endl;

    float sum_indices = 0;
    for (int i = 0; i < singularities_.size(); i++) {
        int vertex_id = singularities_[i].first;
        float index = singularities_[i].second;

        sum_indices += index;

        b(vertex_id) = K(vertex_id) - 2 * M_PI * index;
        std::cout << "Updated defect for vertex " << vertex_id << " from " << K(vertex_id) << " to " << b(vertex_id) << std::endl;

    }

    if ( sum_indices != 2-2*nGenerators ) {
        showMessage("Indices must sum to " + std::to_string(2-2*nGenerators) + " for this mesh!");
        return;
    }

    // std::cout << "Solver: Init." << std::endl;
    // Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
    // Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
    // solver.compute(A);
    // std::cout << "Solver: Computed A." << std::endl;

    // Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >  solver2;
    // solver2.compute(solver.matrixR().transpose());
    // Eigen::MatrixXd y = solver2.solve(solver.colsPermutation().transpose() * (-b));

    // Eigen::MatrixXd x_optimal = solver.matrixQ() * y;


    // One simple solution
    std::cout << "QRSolver: Solving x..." << std::endl;
    Eigen::MatrixXd x = QRSolver.solve(-b);
    std::cout << "QRSolver: Solved x. Start computing d1*d1.T..." << std::endl;

    // Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver2;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver2;
    solver2.compute(d1_ * d1_.transpose());
    std::cout << "solver2: Computed d1*d1.T." << std::endl;

    // Eigen::SparseMatrix<double> I(mesh_->n_faces(),mesh_->n_faces());
    // I.setIdentity();
    Eigen::MatrixXd y = solver2.solve(d1_ * x);
    std::cout << "solver2: Solved (d1*d1.T)^-1. " << std::endl;

    std::cout << "rows " << y.rows() << " cols " << y.cols() << std::endl;

    // Eigen::MatrixXd x_optimal = x - d1_.transpose() * inv * d1_ * x;
    Eigen::MatrixXd x_optimal = x - d1_.transpose() * y;
    std::cout << "Computed x_optimal. " << std::endl;

    for (auto eh : mesh_->edges()) {
        mesh_->property(edge_adjuestment_angle_, mesh_->edge_handle(eh.idx())) = x_optimal(eh.idx());
    }

    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS, DATA_TRIANGLE_MESH); o_it != PluginFunctions::objectsEnd(); ++o_it) {

        auto tri_obj = PluginFunctions::triMeshObject(*o_it);
        emit updatedObject(tri_obj->id(), UPDATE_ALL);
    }
    
    std::cout << "Start building the vector field..." << std::endl;
    showVectorField();

    uint64_t end_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    std::cout << "Completed! Time: " << end_time - start_time << "ms." << std::endl;

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

    // OpenMesh::SmartFaceHandle root_face = *(mesh_->faces_begin());
    OpenMesh::SmartFaceHandle root_face;
    if (constrained_faces_.size() > 0) {
        root_face = constrained_faces_[0];
        mesh_->property(face_tangent_vector_, root_face) = root_face_dir_.normalize();
    } else {
        root_face = *(mesh_->faces_begin());
            OpenMesh::SmartHalfedgeHandle root_heh;

        for (auto heh : root_face.halfedges()) {
            root_heh = heh;
            break;
        }

        ACG::Vec3d init_direction = (mesh_->point(root_heh.to()) - mesh_->point(root_heh.from())).normalize();

        std::cout << "First direction: " << init_direction << " " << root_heh.from().idx() << "->" << root_heh.to().idx() << std::endl;
        auto second_heh = root_heh.next();
        ACG::Vec3d second_direction = (mesh_->point(second_heh.to()) - mesh_->point(second_heh.from())).normalize();
        // std::cout << "Second direction: " << second_direction << " " << second_heh.from().idx() << "->" << second_heh.to().idx() << std::endl;

        // mesh_->property(face_tangent_vector_, root_face) = init_direction.normalize();
        mesh_->property(face_tangent_vector_, root_face) = -( - init_direction + second_direction ).normalize();

    }
    // OpenMesh::SmartFaceHandle root_face = OpenMesh::SmartFaceHandle(420, mesh_);

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
            if (!tri_obj->getAdditionalNode(lineNode, name(), "Vector Field"))
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

                ACG::Vec3d p0 = centroid - scale / 2. * dir;
                ACG::Vec3d p1 = centroid + scale / 2. * dir;

                lineNode->add_color(OpenMesh::Vec4f(0.0f,0.0f,1.0f,1.0f));
                lineNode->add_line(p0, p1);

                lineNode->add_color(OpenMesh::Vec4f(0.0f,1.0f,0.0f,1.0f));
                lineNode->add_line(.2*p0+.8*p1, p1);

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

