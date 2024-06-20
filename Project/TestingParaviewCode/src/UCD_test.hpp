#ifndef __UCD_test_HPP__
#define __UCD_test_HPP__

#include "UCDUtilities.hpp"
#include <gtest/gtest.h>


// ***************************************************************************
TEST(TestUCDUtilities, UCDUtilities_Test0Ds)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points1 = (Eigen::MatrixXd(3, 4)<< 0.8, 1.0, 1.0, 0.8,
                                     0.0, 0.0, 1.0, 1.0,
                                     0.0, 0.0, 0.0, 0.0).finished();
    const Eigen::MatrixXd points2 = (Eigen::MatrixXd(3, 4)<< 0.0, 0.8, 0.8, 0.0,
                                     0.0, 0.0, 0.5, 0.5,
                                     0.0, 0.0, 0.0, 0.0).finished();
    const Eigen::MatrixXd points3 = (Eigen::MatrixXd(3, 4)<< 0.8, 0.8, 0.0, 0.0,
                                     0.5, 1.0, 1.0, 0.5,
                                     0.0, 0.0, 0.0, 0.0).finished();


    exporter.ExportPoints(exportFolder + "/Geometry10Ds.inp",
                          points1);
    exporter.ExportPoints(exportFolder + "/Geometry20Ds.inp",
                          points2);
    exporter.ExportPoints(exportFolder + "/Geometry30Ds.inp",
                          points3);


}


// ***************************************************************************
TEST(TestUCDUtilities, UCDUtilities_Test1Ds)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points1 = (Eigen::MatrixXd(3, 4)<< 0.8, 1.0, 1.0, 0.8,
                                     0.0, 0.0, 1.0, 1.0,
                                     0.0, 0.0, 0.0, 0.0).finished();
    const Eigen::MatrixXd points2 = (Eigen::MatrixXd(3, 4)<< 0.0, 0.8, 0.8, 0.0,
                                     0.0, 0.0, 0.5, 0.5,
                                     0.0, 0.0, 0.0, 0.0).finished();
    const Eigen::MatrixXd points3 = (Eigen::MatrixXd(3, 4)<< 0.8, 0.8, 0.0, 0.0,
                                     0.5, 1.0, 1.0, 0.5,
                                     0.0, 0.0, 0.0, 0.0).finished();
    const Eigen::MatrixXi edges1 = (Eigen::MatrixXi(2, 4)<< 3, 2, 1, 0,
                                    2, 1, 0, 3).finished();
    const Eigen::MatrixXi edges2 = (Eigen::MatrixXi(2, 4)<< 3, 2, 1, 0,
                                    2, 1, 0, 3).finished();
    const Eigen::MatrixXi edges3 = (Eigen::MatrixXi(2, 4)<< 3, 2, 1, 0,
                                    2, 1, 0, 3).finished();

    exporter.ExportSegments(exportFolder + "/Geometry11Ds.inp",
                            points1,
                            edges1);
    exporter.ExportSegments(exportFolder + "/Geometry21Ds.inp",
                            points2,
                            edges2);
    exporter.ExportSegments(exportFolder + "/Geometry31Ds.inp",
                            points3,
                            edges3);
}



/*TEST(TestUCDUtilities, UCDUtilities_Test0Ds)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points1 = (Eigen::MatrixXd(3, 3)<< 0.158305, 0.14992, 0.163228,
                                     0.321895, 0.323377, 0.351296,
                                     0.43477, 0.377919, 0.432974).finished();
    const Eigen::MatrixXd points2 = (Eigen::MatrixXd(3, 4)<< 0.238807, 0.304934, 0.315867, 0.241961,
                                     0.890697, 1.26489, 1.26296, 0.821544,
                                     0.303164, 0.303164, 0.37729, 0.404251).finished();
    const Eigen::MatrixXd points3 = (Eigen::MatrixXd(3, 3)<< 0.254659, 0.210903, 0.251039,
                                     0.543121, 0.451321, 0.622505,
                                     0.811247, 0.63022, 0.695205).finished();
    const Eigen::MatrixXd points4 = (Eigen::MatrixXd(3, 5)<< 0.283902, 0.231674, 0.179673, 0.210903, 0.254659,
                                     0.604473, 0.308929, 0.318119, 0.451321, 0.543121,
                                     0.932231, 0.932231, 0.579654, 0.63022, 0.811247).finished();
    const Eigen::MatrixXd points5 = (Eigen::MatrixXd(3, 6)<< 0.315867, 0.397713, 0.283902, 0.254659, 0.251039, 0.241961,
                                     1.26296, 1.2485, 0.604473, 0.543121, 0.622505, 0.821544,
                                     0.37729, 0.932231, 0.932231, 0.811247, 0.695205, 0.404251).finished();
    const Eigen::MatrixXd points6 = (Eigen::MatrixXd(3, 5)<< 0.138895, 0.204178, 0.196576, 0.163228, 0.14992,
                                     0.325325, 0.694743, 0.550473, 0.351296, 0.323377,
                                     0.303164, 0.303164, 0.420808, 0.432974, 0.377919).finished();
    const Eigen::MatrixXd points7 = (Eigen::MatrixXd(3, 4)<< 0.204178, 0.238807, 0.241961, 0.196576,
                                     0.694743, 0.890697, 0.821544, 0.550473,
                                     0.303164, 0.303164, 0.404251, 0.420808).finished();
    const Eigen::MatrixXd points8 = (Eigen::MatrixXd(3, 5)<< 0.210903, 0.18892, 0.196576, 0.241961, 0.251039,
                                     0.451321, 0.405199, 0.550473, 0.821544, 0.622505,
                                     0.63022, 0.539271, 0.420808, 0.404251, 0.695205).finished();
    const Eigen::MatrixXd points9 = (Eigen::MatrixXd(3, 3)<< 0.18892, 0.163228, 0.196576,
                                     0.405199, 0.351296, 0.550473,
                                     0.539271, 0.432974, 0.420808).finished();
    const Eigen::MatrixXd points10 = (Eigen::MatrixXd(3, 5)<< 0.179673, 0.158305, 0.163228, 0.18892, 0.210903,
                                     0.318119, 0.321895, 0.351296, 0.405199, 0.451321,
                                     0.579654, 0.43477, 0.432974,  0.539271, 0.63022).finished();

    exporter.ExportPoints(exportFolder + "/Geometry10Ds.inp",
                          points1);
    exporter.ExportPoints(exportFolder + "/Geometry20Ds.inp",
                          points2);
    exporter.ExportPoints(exportFolder + "/Geometry30Ds.inp",
                          points3);
    exporter.ExportPoints(exportFolder + "/Geometry40Ds.inp",
                          points4);
    exporter.ExportPoints(exportFolder + "/Geometry50Ds.inp",
                          points5);
    exporter.ExportPoints(exportFolder + "/Geometry60Ds.inp",
                          points6);
    exporter.ExportPoints(exportFolder + "/Geometry70Ds.inp",
                          points7);
    exporter.ExportPoints(exportFolder + "/Geometry80Ds.inp",
                          points8);
    exporter.ExportPoints(exportFolder + "/Geometry90Ds.inp",
                          points9);
    exporter.ExportPoints(exportFolder + "/Geometry100Ds.inp",
                          points10);

}

// ***************************************************************************

TEST(TestUCDUtilities, UCDUtilities_Test1Ds)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points1 = (Eigen::MatrixXd(3, 3)<< 0.158305, 0.14992, 0.163228,
                                     0.321895, 0.323377, 0.351296,
                                     0.43477, 0.377919, 0.432974).finished();
    const Eigen::MatrixXd points2 = (Eigen::MatrixXd(3, 4)<< 0.238807, 0.304934, 0.315867, 0.241961,
                                     0.890697, 1.26489, 1.26296, 0.821544,
                                     0.303164, 0.303164, 0.37729, 0.404251).finished();
    const Eigen::MatrixXd points3 = (Eigen::MatrixXd(3, 3)<< 0.254659, 0.210903, 0.251039,
                                     0.543121, 0.451321, 0.622505,
                                     0.811247, 0.63022, 0.695205).finished();
    const Eigen::MatrixXd points4 = (Eigen::MatrixXd(3, 5)<< 0.283902, 0.231674, 0.179673, 0.210903, 0.254659,
                                     0.604473, 0.308929, 0.318119, 0.451321, 0.543121,
                                     0.932231, 0.932231, 0.579654, 0.63022, 0.811247).finished();
    const Eigen::MatrixXd points5 = (Eigen::MatrixXd(3, 6)<< 0.315867, 0.397713, 0.283902, 0.254659, 0.251039, 0.241961,
                                     1.26296, 1.2485, 0.604473, 0.543121, 0.622505, 0.821544,
                                     0.37729, 0.932231, 0.932231, 0.811247, 0.695205, 0.404251).finished();
    const Eigen::MatrixXd points6 = (Eigen::MatrixXd(3, 5)<< 0.138895, 0.204178, 0.196576, 0.163228, 0.14992,
                                     0.325325, 0.694743, 0.550473, 0.351296, 0.323377,
                                     0.303164, 0.303164, 0.420808, 0.432974, 0.377919).finished();
    const Eigen::MatrixXd points7 = (Eigen::MatrixXd(3, 4)<< 0.204178, 0.238807, 0.241961, 0.196576,
                                     0.694743, 0.890697, 0.821544, 0.550473,
                                     0.303164, 0.303164, 0.404251, 0.420808).finished();
    const Eigen::MatrixXd points8 = (Eigen::MatrixXd(3, 5)<< 0.210903, 0.18892, 0.196576, 0.241961, 0.251039,
                                     0.451321, 0.405199, 0.550473, 0.821544, 0.622505,
                                     0.63022, 0.539271, 0.420808, 0.404251, 0.695205).finished();
    const Eigen::MatrixXd points9 = (Eigen::MatrixXd(3, 3)<< 0.18892, 0.163228, 0.196576,
                                     0.405199, 0.351296, 0.550473,
                                     0.539271, 0.432974, 0.420808).finished();
    const Eigen::MatrixXd points10 = (Eigen::MatrixXd(3, 5)<< 0.179673, 0.158305, 0.163228, 0.18892, 0.210903,
                                      0.318119, 0.321895, 0.351296, 0.405199, 0.451321,
                                      0.579654, 0.43477, 0.432974,  0.539271, 0.63022).finished();

    const Eigen::MatrixXi edges1 = (Eigen::MatrixXi(2, 3)<< 2, 1, 0,
                                    1, 0, 2).finished();
    const Eigen::MatrixXi edges2 = (Eigen::MatrixXi(2, 4)<< 3, 2, 1, 0,
                                    2, 1, 0, 3).finished();
    const Eigen::MatrixXi edges3 = (Eigen::MatrixXi(2, 3)<< 2, 1, 0,
                                    1, 0, 2).finished();
    const Eigen::MatrixXi edges4 = (Eigen::MatrixXi(2, 5)<< 4, 3, 2, 1, 0,
                                    3, 2, 1, 0, 4).finished();
    const Eigen::MatrixXi edges5 = (Eigen::MatrixXi(2, 6)<< 5, 4, 3, 2, 1, 0,
                                    4, 3, 2, 1, 0, 5).finished();
    const Eigen::MatrixXi edges6 = (Eigen::MatrixXi(2, 5)<< 4, 3, 2, 1, 0,
                                    3, 2, 1, 0, 4).finished();
    const Eigen::MatrixXi edges7 = (Eigen::MatrixXi(2, 4)<< 3, 2, 1, 0,
                                    2, 1, 0, 3).finished();
    const Eigen::MatrixXi edges8 = (Eigen::MatrixXi(2, 5)<< 4, 3, 2, 1, 0,
                                    3, 2, 1, 0, 4).finished();
    const Eigen::MatrixXi edges9 = (Eigen::MatrixXi(2, 3)<< 2, 1, 0,
                                    1, 0, 2).finished();
    const Eigen::MatrixXi edges10 = (Eigen::MatrixXi(2, 5)<< 4, 3, 2, 1, 0,
                                    3, 2, 1, 0, 4).finished();

    exporter.ExportSegments(exportFolder + "/Geometry11Ds.inp",
                            points1,
                            edges1);
    exporter.ExportSegments(exportFolder + "/Geometry21Ds.inp",
                            points2,
                            edges2);
    exporter.ExportSegments(exportFolder + "/Geometry31Ds.inp",
                            points3,
                            edges3);
    exporter.ExportSegments(exportFolder + "/Geometry41Ds.inp",
                            points4,
                            edges4);
    exporter.ExportSegments(exportFolder + "/Geometry51Ds.inp",
                            points5,
                            edges5);
    exporter.ExportSegments(exportFolder + "/Geometry61Ds.inp",
                            points6,
                            edges6);
    exporter.ExportSegments(exportFolder + "/Geometry71Ds.inp",
                            points7,
                            edges7);
    exporter.ExportSegments(exportFolder + "/Geometry81Ds.inp",
                            points8,
                            edges8);
    exporter.ExportSegments(exportFolder + "/Geometry91Ds.inp",
                            points9,
                            edges9);
    exporter.ExportSegments(exportFolder + "/Geometry101Ds.inp",
                            points10,
                            edges10);

}*/

// ***************************************************************************
TEST(TestUCDUtilities, UCDUtilities_Test2Ds)
{
    std::string exportFolder = "./";

    Gedim::UCDUtilities exporter;
    const Eigen::MatrixXd points = (Eigen::MatrixXd(3, 8)<< 0.8, 1.0, 1.0, 0.8, 0.0, 0.8, 0.0, 0.0,
                                    0.0, 0.0, 1.0, 1.0, 0.0, 0.5, 0.5, 1.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0).finished();
    const std::vector<std::vector<unsigned int>> polygons =
        {
            { 0, 1, 2, 3 },
            { 4, 0, 5, 6 },
            { 5, 3, 7, 6 }
        };

    exporter.ExportPolygons(exportFolder + "/Geometry2Ds.inp",
                            points,
                            polygons);
}
// ***************************************************************************

#endif // __UCD_test_HPP__
