#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QVector3D>
#include <QLabel>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

/**
  * @brief Données de résultat de la persistance
  */
typedef struct {
    int dim;                /// dimension of hole
    float t_radius;         /// radius of the thickness ball
    float t_x, t_y, t_z;    /// center of the T ball
    float b_radius;         /// radius of the breadth ball
    float b_x, b_y, b_z;    /// center of the B ball
    float persistence;    /// t_radius + b_radius
    bool present_hole;
    //bool infinite_t;
    //bool infinite_b;
} TBball;

/**
 * @brief Classe principale de l'application
 */
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    /**
     * @brief Création de la fenêtre
     */
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1) const;
//    void displayPath(std::vector<CPoint3D> resultpoints);
    void resetAllColorsAndThickness(MyMesh* _mesh);

    void print_tb_pairs();

private slots:
    void on_pushButton_test_clicked();
    void on_checkBox_wireframe_stateChanged();
    void on_checkBox_allballs_stateChanged();
    void on_checkBox_presentholes_stateChanged();
    void on_comboBox_TBchoice_currentIndexChanged();
    void on_spinBox_hole_valueChanged(int i);
    void on_checkBox_dim0_stateChanged();
    void on_checkBox_dim1_stateChanged();
    void on_checkBox_dim2_stateChanged();

    void on_slider_min_persistence_valueChanged();
    void on_spinbox_min_persistence_valueChanged();
    //void update_min_persistence(double value);
public slots:
    void read_mesh();
    void read_tb();
    void selectPair(int dim, int pos);


private:

    /**
     * @brief Maillage utilisant la structure d'OpenMesh
     */
    MyMesh mesh;
    std::vector<TBball> m_tb_balls;
    std::vector<bool> m_displayed_balls;
    bool m_thickness; // show thickness ball (breadth, otherwise)
    float max_persistence;
    float current_min_persistence;

    QString m_filebasename; // basename of the filename of the mesh

    Ui::MainWindow *ui;

/* Methods */
    void display_balls();
    void show_tb_balls(const std::vector<TBball>& balls);
};


#endif // MAINWINDOW_H
