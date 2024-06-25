#ifndef BUNCH_H
#define BUNCH_H

#include <algorithm>
#include <exception>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>


#include "wrap.h"

namespace ris2{

/// @brief Classes for working with the number of drawable objects all at ones.

/// @brief A class for handling the painting job
class Palette{
public:
  /// @brief An exception thrown in the case of number of number of objects to paint is greater then the number of styles availiable
  class VectorOfStylesIsOutOfRange : public std::exception{
  public:
    VectorOfStylesIsOutOfRange( size_t req_idx) {
      error_message_ = std::string("Required index ").append( std::to_string(req_idx) ).append(" is out of range of vector of styles");
    }
    const char* what() const noexcept override {
      return error_message_.c_str();
    }
  private:
    std::string error_message_{};
  };
  Palette() = default;
  ~Palette() = default;
  
  Palette& SetPalette( std::vector<Style> palette ){ styles_ = std::move(palette); return *this; }
  const std::vector<Style>& GetPalette(){ return styles_; }
  template< class T >
  void PaintObjects( std::vector< Wrap<T> >& objects ) const { 
    std::for_each( objects.begin(), objects.end(), 
    [this, i=0]( auto& obj ) mutable {
      try{
        obj.SetStyle( styles_.at(i) );
      } catch( std::out_of_range ){
        throw VectorOfStylesIsOutOfRange(i);
      }
      ++i;
    } ); 
  }
  void PaintObject( size_t style_idx, TGraphErrors* graph ) { 
    styles_.at(style_idx).operator()( graph );
  }
  TLegend* MakeLegend( const std::vector<std::string>& captions, std::vector<double> position = {} ){
    int i=0;
    TLegend* leg{};
    if( position.size() == 4 )
      leg = new TLegend{ position.at(0), position.at(1), position.at(2), position.at(3) };
    else{
      leg = new TLegend{};
    }
    for( auto& titile : captions ){
      auto graph = new TGraph(0);
      graph->SetLineColor( styles_.at(i).color_ );
      graph->SetMarkerColor( styles_.at(i).color_ );
      if( styles_.at(i).marker_ >= 0 ){
        graph->SetMarkerStyle( styles_.at(i).marker_ );
        leg->AddEntry( graph, titile.c_str(), "P" );
      }
      if( styles_.at(i).marker_ < 0 ){
        graph->SetLineStyle( styles_.at(i).marker_ );
        leg->AddEntry( graph, titile.c_str(), "L" );
      }
      ++i;
    }
    return leg;
  }
protected:
  std::vector<Style> styles_{};
};

/// @brief A class wrapping the std::vector<Wrap<T>>.
template<typename T>
class Bunch {
public:
  /// @brief Default constructor
  Bunch<T>() = default;
  /// @brief Copy constructor
  Bunch<T>( const Bunch<T>& other) = default;
  /// @brief Move constructor
  Bunch<T>( Bunch<T>&& other) = default;
  /// @brief Copy assignment operator
  Bunch<T>& operator=( const Bunch<T>& other) = default;
  /// @brief Move assignment operator
  Bunch<T>& operator=( Bunch<T>&& ) = default;
  /// @brief Default destructor
  ~Bunch<T>() = default;

  /// @brief Constructs and places the value to the end of the underlying vector
  /// @param Args are forwarded to the constructor of Wrap<T>. See the documentation for Wrap<T>
  template<typename... Args>
  Bunch& AddToBunch( std::string title, Args... args ){ 
    bunch_.emplace_back( title, args... ); 
    return *this;
  }
  /// @brief Getter for the palette
  Palette& GetPalette() { return palette_; }
  /// @brief Function to perform an action on all the ellements of the underlyind std::vector<Wrap<T>>
  /// @param function function to perform on each Wrap<T>. The function should have void return type
  template<typename Func>
  Bunch& Perform( Func function ){ std::for_each( bunch_.begin(), bunch_.end(), function ); return *this; }
  /// @brief Function to access the underlying vector of Wrap<T>
  std::vector<Wrap<T>>& operator*(){
    palette_.PaintObjects(bunch_);
    return bunch_; 
  }
  /// @brief Same as Perform. Added for brievity. 
  /// @param function function to perform on each Wrap<T>. The function should have void return type
  template<typename Func>
  Bunch& operator()( Func function ){ std::for_each( bunch_.begin(), bunch_.end(), function ); return *this; }
  /// @brief Function for index access to the underlying vector of Wrap<T>
  Wrap<T>& operator[]( size_t idx ){ return bunch_.at(idx); }
  /// @brief Creates the legend using the titles and styles of the Wrap<T> placed to the vector.
  /// @param position Coordinates of the bottom-left and up-right corners of the legend on the canvas { x1, y1, x2, y2 }. 
  /// Relative positioning is used (i.e. coordinates from 0 to 1). 
  TLegend* MakeLegend(std::vector<double> position = {}){
    std::vector<std::string> bunch_titles_;
    std::for_each( bunch_.begin(), bunch_.end(), [&bunch_titles_]( const Wrap<T>& obj ) mutable { 
      bunch_titles_.emplace_back( obj.GetTitle() ); 
    } );
    auto leg = palette_.MakeLegend( bunch_titles_, position );
    return leg;
  }

private:
  Palette palette_{};
  std::vector< Wrap<T> > bunch_{};
};

/// @brief A class to automatically plot double-differential results
/// Draws the projectin on the required axis in narrow bins of the second axis

#ifdef USE_QNTOOLS

class DoubleDifferential {
public:
  /// @brief Default constructor
  DoubleDifferential() = delete;
  /// @brief Constructor initializing the Wrap<Correlation> which will be sliced
  /// @param title can be leaved empty. 
  /// @param Args is a param pack for Result<Correlation> construction
  template<typename ...Args>
  DoubleDifferential(std::string title, Args... args) : base_correlation_{ title, args... }{ }
  /// @brief Copy constructor
  DoubleDifferential( DoubleDifferential& ) = default;
  /// @brief Copy assignment operator
  DoubleDifferential& operator=( DoubleDifferential& ) = default;
  /// @brief Move constructor
  DoubleDifferential( DoubleDifferential&& ) = default;
  /// @brief Move assignment operator
  DoubleDifferential& operator=( DoubleDifferential&& ) = default;
  /// @brief Default destructor
  ~DoubleDifferential() = default;
  /// @brief Function to perform on the underlying Wrap<Correlation>
  /// @param func should have the return-type void
  template<typename Func>
  DoubleDifferential& Perform( const Func& function ){ 
    function( base_correlation_ ); 
    return *this; 
  }
  /// @brief Same function as perform for brievity 
  template<typename Func>
  DoubleDifferential& operator()( const Func& function ){ 
    function( base_correlation_ ); 
    return *this; 
  }
  Palette& GetPalette(){ return palette_; }
  DoubleDifferential& SetPalette( const std::vector<Style>& styles ){ palette_.SetPalette(styles); return *this; }
  /// @brief Setter for the axis in bins of which the correlation will be projected on the projection axis 
  DoubleDifferential& SetSliceAxis(Qn::AxisD slie_axis){ slice_axis_ = std::move(slie_axis); return *this; }
  /// @brief Setter for the axis on which the correlation will be projected in the bins of slice axis 
  DoubleDifferential& SetProjectionAxis(Qn::AxisD projection_axis){ projection_axis_ = std::move(projection_axis); return *this; }
  /// @brief Creates the legend out of all the slices.
  /// @param slice_var text description of the values on the slice axis
  /// @param slice_var units of measure of the slice axis
  /// @param position Coordinates of the bottom-left and up-right corners of the legend on the canvas { x1, y1, x2, y2 }. 
  /// Relative positioning is used (i.e. coordinates from 0 to 1). 
  TLegend* MakeLegend(const std::string& slice_var, const std::string& units="", std::vector<double> position = {} ){
    auto bin_edges = slice_axis_.GetBinEdges();
    std::vector<std::string> captions{};
    for( int i = 0; i < bin_edges.size() - 1; ++i ){
      auto stream = std::ostringstream{};
      stream << std::setprecision(2) << slice_var << " " << bin_edges.at(i) << "-" << bin_edges.at(i+1) << " " << units;
      captions.emplace_back( stream.str() );
    }
    return palette_.MakeLegend( captions, position );
  }
  /// @brief returns the result slices
  std::vector<Wrap<Correlation>>& operator*(){
    if( result_correlations_.empty() )
      Slice();
    return result_correlations_;
  }
private:
  void Slice(){
    auto n_projections = slice_axis_.GetNBins();
    auto slice_var = slice_axis_.Name();
    auto bin_edges = slice_axis_.GetBinEdges();
    // auto projection_var = projection_axis_.Name();
    base_correlation_
      .Rebin( std::vector<Qn::AxisD>{slice_axis_, projection_axis_} )
      .Project( std::vector<Qn::AxisD>{slice_axis_, projection_axis_} );
    result_correlations_.reserve(n_projections);
    for( int i=0; i<n_projections; ++i ){
      result_correlations_.emplace_back( base_correlation_ );
      
      result_correlations_.back().Rebin( std::vector<Qn::AxisD>{{slice_var, 1, bin_edges.at(i), bin_edges.at(i+1)} } )
        .Project( std::vector<Qn::AxisD>{{projection_axis_}} );
      palette_.PaintObjects( result_correlations_ );
    }
  }
  Palette palette_;
  Wrap<Correlation> base_correlation_;
  std::vector< Wrap<Correlation> > result_correlations_;
  Qn::AxisD slice_axis_{};
  Qn::AxisD projection_axis_{};
};

#endif

/// @brief Class to create the ratio plots.
/// Used for systematic checks. The ratio is built with respect to one of the selected values.
template<typename T>
class RatioBuilder {
public:
  /// @brief default constructor
  RatioBuilder<T>() = delete;
  /// @brief Constructor intializing the reference value with respect to which the ratios will be built.
  /// @param title is the title used for creating the legend
  /// @param args is a parameter pack to initialize the underlying Result<T>
  template<typename... Args>
  RatioBuilder( std::string title, Args... args ) {
    results_.emplace_back( args... );
    titles_.emplace_back( std::move(title) );
  }
  /// @brief Default constuctor
  ~RatioBuilder() = default;
  /// @brief Function to add the object which the ratio will be built for
  /// @param title is the title used for creating the legend
  /// @param args is the parameter pack to initialize the Result<T>
  template<typename ...Args>
  RatioBuilder& AddToBunch( std::string title, Args... args ){ 
    titles_.emplace_back( std::move(title) );
    results_.emplace_back( args...); 
    return *this; 
  }
  /// @brief Returns the Result<T> objects for which the ratios are built
  std::vector<Result<T>>& GetResults(){
    return results_;
  }
  /// @brief Returns the Result<T> objects containing the ratios
  std::vector<Result<T>>& GetRatios(){
    BuildRatios();
    return ratios_;
  }
  /// @brief Function to perform the action (projection/rebin/scale, etc.) on all the objects for which the ratios will be built.
  /// @param function accepts Ratio<T>& as a parameter and returns nothing.
  template<typename Func>
  RatioBuilder<T>& Perform( const Func& function ){ 
    std::for_each( results_.begin(), results_.end(), function ); 
    return *this; 
  }
  /// @brief Same as Perform-function.
  template<typename Func>
  RatioBuilder<T>& operator()( const Func& function ){ 
    std::for_each( results_.begin(), results_.end(), function ); 
    return *this; 
  }
  RatioBuilder<T>& SetPalette( std::vector<Style> corr_styles ){
    palette_.SetPalette(corr_styles);
    return *this;
  }
  /// @brief Returns the reference Result<T> relative to which the ratios are built 
  Result<T>& GetReference(){ return results_.front(); }
  /// @brief Creates the legend out of all the slices.
  /// @param position Coordinates of the bottom-left and up-right corners of the legend on the canvas { x1, y1, x2, y2 }. 
  /// Relative positioning is used (i.e. coordinates from 0 to 1). 
  TLegend* MakeLegend( std::vector<double> position = {} ){
    auto leg = palette_.MakeLegend(titles_, position);
    return leg;
  }
  /// @brief Returns the number of Ratio<T> for which the ratios will be built 
  size_t Size() const { return results_.size(); }
  /// @brief Returns the points for which the ratios are built 
  std::vector<Wrap<TGraphErrors>>& GetResultWraps(){ UpdatePoints(); return result_points_; }
  /// @brief Returns the points of the ratios 
  std::vector<Wrap<TGraphErrors>>& GetRatioWraps(){ UpdatePoints(); return ratio_points_; }
private:
  void BuildRatios(){
    if( ratios_.empty() ){
      int i=0;
      const auto& reference = results_.front();
      ratios_.reserve( results_.size() );
      std::for_each( results_.begin(), results_.end(), [&reference, this, i=0](auto& obj) mutable {
        ratios_.emplace_back( obj.Divide( reference )  );
      } );
    }
  }
  void UpdatePoints(){
    BuildRatios();
    if( ! result_points_.empty() ){
      return;
    }
    std::for_each( results_.begin(), results_.end(), 
    [this]( auto obj  ) {
      result_points_.push_back( Wrap<TGraphErrors>(obj.GetPoints()) );
    });
    std::for_each( ratios_.begin(), ratios_.end(), 
    [this]( auto obj  ) {
      ratio_points_.push_back( Wrap<TGraphErrors>(obj.GetPoints()) );
    });
    palette_.PaintObjects(result_points_);
    palette_.PaintObjects(ratio_points_);
  };
  std::vector<std::string>  titles_{};
  Palette palette_;
  std::vector<Result<T>> results_{};
  std::vector<Result<T>> ratios_{};
  std::vector<Wrap<TGraphErrors>> result_points_{};
  std::vector<Wrap<TGraphErrors>> ratio_points_{};
};

}

#endif // BUNCH_H
