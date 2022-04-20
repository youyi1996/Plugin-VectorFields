#include "vectorFields.hh"
#include "OpenFlipper/BasePlugin/PluginFunctions.hh"

#include <string>

void VectorFields::initializePlugin() {
    QWidget* toolBox = new QWidget();

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

    QGridLayout* layout = new QGridLayout(toolBox);
    layout->addWidget(singularityLabel, 0, 0, 1, 3);
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


    // connect(smoothButton, SIGNAL(clicked()), this, SLOT(simpleLaplace()));

    connect(addSingularityButton, SIGNAL(clicked()), this, SLOT(addSingularity()));
    connect(removeSingularityButton, SIGNAL(clicked()), this, SLOT(removeSingularity()));
    connect(clearSingularityButton, SIGNAL(clicked()), this, SLOT(clearSingularities()));

    connect(singularityTable->selectionModel(), SIGNAL(selectionChanged(const QItemSelection&, const QItemSelection&)), this, SLOT(onItemSelected(const QItemSelection&, const QItemSelection&)));

    connect(singularityTable, SIGNAL(cellChanged(int, int)), this, SLOT(onCellChanged(int, int)));


    connect(runButton, SIGNAL(clicked()), this, SLOT(runAll()));

    emit addToolbox(tr("VectorField"), toolBox);

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

            mesh_ = trimesh;

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

    
}