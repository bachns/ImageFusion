/********************************************************************************
*   Copyright (C) 2018 by Bach Nguyen Sy                                       *
*   Unauthorized copying of this file via any medium is strictly prohibited    *
*                                                                              *
*   <bachns.dev@gmail.com>                                                     *
*   https://bachns.wordpress.com                                               *
*   https://www.facebook.com/bachns.dev                                        *
*                                                                              *
********************************************************************************/

/**
* File name:    Fusion/Fusion.cpp
* Date created: Monday, Apr 22, 2019
* Written by Bach Nguyen Sy
*/

#include "otbVectorImage.h"
#include "otbImageFileReader.h"
#include "otbImageFileWriter.h"
#include "otbStreamingResampleImageFilter.h"
#include "otbGridResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "otbBCOInterpolateImageFunction.h"
#include "otbStandardFilterWatcher.h"
#include "otbSimpleRcsPanSharpeningFusionImageFilter.h"
#include "boost/program_options.hpp"


int main(int argc, char* argv[]) {
	namespace po = boost::program_options;
	po::options_description desc("The Image Fusion Program written by Bach NS\nAllowed options");
	desc.add_options()
		("help,h", "Show help")
		("pan-image,p", po::value<std::string>(), "Set the panchromatic image")
		("xs-image,x", po::value<std::string>(), "Set the multispectral image")
		("output-image,o", po::value<std::string>()->default_value("a.tif"), "Set the output image\n  [default=a.tif]")
		("interpolator,i", po::value<std::string>()->default_value("ln"), "Set the interpolation method:\n  * nn = Nearest neighbor\n  * ln = Linear\n  * bc = Bicubic\n  [default=ln]")
		("radius,r", po::value<int>()->default_value(2), "Set radius for bicubic interpolation\n  [default=2]");
	po::variables_map vm;

	std::string pan_file_name;
	std::string xs_file_name;
	std::string output_file_name;
	std::string interp_method;
	int radius;

	try {
		store(parse_command_line(argc, argv, desc), vm);
		notify(vm);
		if (vm.count("help")) {
			std::cout << desc << std::endl;
			return EXIT_SUCCESS;
		}
		if (!vm.count("pan-image") || !vm.count("xs-image")) {
			throw po::error("The option '--pan-image' or '--xs-image' is required but missing");
		}
		pan_file_name = vm["pan-image"].as<std::string>();
		xs_file_name = vm["xs-image"].as<std::string>();
		output_file_name = vm["output-image"].as<std::string>();
		interp_method = vm["interpolator"].as<std::string>();
		radius = vm["radius"].as<int>();
	}
	catch (po::error& e) {
		std::cerr << "Error: " << e.what() << std::endl << std::endl;
		std::cerr << desc << std::endl;
		return EXIT_FAILURE;
	}

	const int Dimension = 2;
	typedef float PixelType;
	typedef otb::Image<PixelType, Dimension> FloatImageType;
	typedef otb::VectorImage<PixelType, Dimension> FloatVectorImageType;
	typedef otb::ImageFileReader<FloatImageType> ImageReaderType;
	typedef otb::ImageFileReader<FloatVectorImageType> VectorImageReaderType;
	typedef otb::ImageFileWriter<FloatVectorImageType> VectorImageWriterType;
	typedef otb::GridResampleImageFilter<FloatVectorImageType, FloatVectorImageType> GridResampleFilterType;


	ImageReaderType::Pointer pan_reader = ImageReaderType::New();
	pan_reader->SetFileName(pan_file_name);
	pan_reader->Update();
	FloatImageType::Pointer pan_image = pan_reader->GetOutput();

	VectorImageReaderType::Pointer xs_reader = VectorImageReaderType::New();
	xs_reader->SetFileName(xs_file_name);
	xs_reader->Update();
	FloatVectorImageType::Pointer xs_image = xs_reader->GetOutput();

	FloatImageType::SpacingType pan_spacing = pan_image->GetSpacing();
	FloatVectorImageType::SpacingType xs_spacing = xs_image->GetSpacing();
	double scale[2];
	scale[0] = pan_spacing[0] / xs_spacing[0];
	scale[1] = pan_spacing[1] / xs_spacing[1];

	if (scale[0] > 1.0 || scale[1] > 1.0) {
		std::cerr << "The multispectral image is better than panchromatic image" << std::endl;
		return EXIT_FAILURE;
	}

	GridResampleFilterType::Pointer grid_resampler = GridResampleFilterType::New();
	grid_resampler->SetInput(xs_image);
	grid_resampler->SetOutputParametersFromImage(xs_image);

	FloatVectorImageType::SpacingType resp_spacing;
	resp_spacing[0] = xs_spacing[0] * scale[0];
	resp_spacing[1] = xs_spacing[1] * scale[1];
	grid_resampler->SetOutputSpacing(resp_spacing);

	FloatVectorImageType::PointType xs_origin = xs_image->GetOrigin();
	FloatVectorImageType::PointType resp_origin;
	resp_origin[0] = xs_origin[0] + 0.5 * xs_spacing[0] * (scale[0] - 1.0);
	resp_origin[1] = xs_origin[1] + 0.5 * xs_spacing[1] * (scale[1] - 1.0);
	grid_resampler->SetOutputOrigin(resp_origin);

	if (!interp_method.compare("nn")) {
		typedef itk::NearestNeighborInterpolateImageFunction<FloatVectorImageType> NearestNeighborInterpolationType;
		NearestNeighborInterpolationType::Pointer interpolator = NearestNeighborInterpolationType::New();
		grid_resampler->SetInterpolator(interpolator);
	}
	else if (!interp_method.compare("bc")) {
		typedef otb::BCOInterpolateImageFunction<FloatVectorImageType> BCOInterpolationType;
		BCOInterpolationType::Pointer interpolator = BCOInterpolationType::New();
		interpolator->SetRadius(radius);
		grid_resampler->SetInterpolator(interpolator);
	}
	else {
		typedef itk::LinearInterpolateImageFunction<FloatVectorImageType> LinearInterpolationType;
		LinearInterpolationType::Pointer interpolator = LinearInterpolationType::New();
		grid_resampler->SetInterpolator(interpolator);
	}

	GridResampleFilterType::SizeType recomputed_size;
	recomputed_size[0] = xs_image->GetLargestPossibleRegion().GetSize()[0] / scale[0];
	recomputed_size[1] = xs_image->GetLargestPossibleRegion().GetSize()[1] / scale[1];
	grid_resampler->SetOutputSize(recomputed_size);

	FloatVectorImageType::PixelType default_value;
	itk::NumericTraits<FloatVectorImageType::PixelType>::SetLength(default_value, xs_image->GetNumberOfComponentsPerPixel());
	grid_resampler->SetEdgePaddingValue(default_value);

	FloatVectorImageType::Pointer resp_image = grid_resampler->GetOutput();

	typedef otb::SimpleRcsPanSharpeningFusionImageFilter<FloatImageType, FloatVectorImageType, FloatVectorImageType> FusionFilterType;
	FusionFilterType::Pointer fusion = FusionFilterType::New();
	fusion->SetPanInput(pan_image);
	fusion->SetXsInput(resp_image);

	VectorImageWriterType::Pointer image_writer = VectorImageWriterType::New();
	image_writer->SetFileName(output_file_name);
	image_writer->SetInput(fusion->GetOutput());
	image_writer->SetAutomaticTiledStreaming();
	otb::StandardFilterWatcher watcher(image_writer, "Image Saving");
	image_writer->Update();
	return EXIT_SUCCESS;
}