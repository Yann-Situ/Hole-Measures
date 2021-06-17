#ifndef CHARTVIEW_H
#define CHARTVIEW_H

#include <QtCharts/QChartView>
#include <QtCharts/QScatterSeries>

QT_CHARTS_USE_NAMESPACE

class ChartView : public QChartView
{
    Q_OBJECT
public:
    explicit ChartView(QWidget *parent = 0);
    void load_pairs(std::vector<int> dim, std::vector<double> x, std::vector<double> y);
    qreal dist(const QPointF &p1, const QPointF &p2) const;

signals:
    void pairSelected(int dim, int pos);

protected:
    void keyPressEvent(QKeyEvent *event);

private Q_SLOTS:
    void handleClickedPoint(const QPointF &point);

private:
    QScatterSeries *m_pairs_0;
    QScatterSeries *m_pairs_1;
    QScatterSeries *m_pairs_2;
};

#endif // CHARTVIEW_H
