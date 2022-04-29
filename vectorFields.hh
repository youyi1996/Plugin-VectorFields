#ifndef VECTORFIELDS_HH_INCLUDED
#define VECTORFIELDS_HH_INCLUDED
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/LoggingInterface.hh>
#include <OpenFlipper/BasePlugin/LoadSaveInterface.hh>
#include <OpenFlipper/BasePlugin/MouseInterface.hh>
#include <OpenFlipper/BasePlugin/PickingInterface.hh>
#include <OpenFlipper/common/Types.hh>
#include <ACG/Scenegraph/LineNode.hh>

#include <QWidget>
#include <QPushButton>
#include <QLabel>
#include <QGridLayout>
#include <QSpinBox>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QHeaderView>
#include <QItemSelectionModel>
#include <QLineEdit>

#include <vector>
#include <map>


#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include "./eigen3/Eigen/Sparse"


class VectorFields : public QObject, BaseInterface, ToolboxInterface, LoggingInterface, LoadSaveInterface, MouseInterface, PickingInterface {
Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)
    Q_INTERFACES(LoadSaveInterface)
    Q_INTERFACES(MouseInterface)
    Q_INTERFACES(PickingInterface)
    Q_PLUGIN_METADATA(IID "org.OpenFlipper.Plugins.examples.VectorFields")


signals:
    void updateView();

    //LoggingInterface
    void log(Logtype _type, QString _message);
    void log(QString _message);

    // ToolboxInterface
    void addToolbox(QString _name, QWidget* _widget);

    // LoadSaveInterface
    void addEmptyObject(DataType _type, int& _id);
    void updatedObject(int _id, const UpdateType& _type);

    // PickingInterface
    void addPickMode(const std::string &_mode);

public:

    typedef std::vector<OpenMesh::HalfedgeHandle> Cycle;
    
    ~VectorFields() {};
    
    QString name() { return QString("VectorFields"); };
    
    QString description() { return QString("Generates a tangent vector field with trivial connections."); };

private:

    QTableWidget* singularityTable;

    TriMesh* mesh_ = nullptr;
    
    bool enableOnCellChanged = true;

    OpenMesh::EPropHandleT<double> edge_adjuestment_angle_;
    OpenMesh::FPropHandleT<ACG::Vec3d> face_tangent_vector_;
    OpenMesh::VPropHandleT<OpenMesh::VertexHandle> parent_prime_;
    OpenMesh::FPropHandleT<OpenMesh::FaceHandle> parent_dual_;

    // A vector that stores singular vertices with indices. <vertex_id, index>.
    std::vector<std::pair<int, float> > singularities_; 

private slots:
    // BaseInterface
    ACG::Vec3d transport(ACG::Vec3d, OpenMesh::SmartFaceHandle&, OpenMesh::SmartFaceHandle&);

    void initializePlugin();

    void updateVertexColors(TriMesh &);
    void initializeMesh();

    void addSingularity();
    void removeSingularity();
    void clearSingularities();

    void refreshSingularityTable();
    void onItemSelected(const QItemSelection&, const QItemSelection&);
    void onCellChanged(int, int);
    

    // Jordan's Part
    bool inPrimalSpanningTree(TriMesh& mesh_, OpenMesh::HalfedgeHandle he);
    bool inDualSpanningTree(TriMesh& mesh_, OpenMesh::HalfedgeHandle he);
    void buildPrimalSpanningTree(TriMesh& mesh_);
    void buildDualSpanningCoTree(TriMesh& mesh_);
    void buildTreeCotreeDecomposition(TriMesh& mesh_);

    OpenMesh::HalfedgeHandle sharedHalfEdge(TriMesh& mesh_, OpenMesh::VertexHandle v, OpenMesh::VertexHandle w);
    OpenMesh::HalfedgeHandle sharedHalfEdge(TriMesh& mesh_, OpenMesh::FaceHandle f, OpenMesh::FaceHandle g);
    bool isDualBoundaryLoop(TriMesh& mesh_, const Cycle& cycle);
    void appendDualGenerators(TriMesh& mesh_, std::vector<Cycle>& cycles);
    // Jordan's Part... END


    void findContractibleLoops(TriMesh& _mesh, std::vector<std::vector<int> >);
    void findNonContractibleLoops(TriMesh& _mesh, std::vector<std::vector<int> >);

    void runAll();

    void showPrimalTree();
    void showDualTree();
    void showCycles(std::vector<Cycle> &);
    void showVectorField();
    void generateFakeVectorField();

public slots:
    QString version() { return QString("1.0"); };
};
#endif