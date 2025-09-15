from .views.variant_view import VariantMetricsView, VariantAnnotationView
from .views.overview_view import OverviewDataView, TopHitsView
from .views.trait_view import ManhattanView, QQView, TraitView, ChromosomeBoundsView, TraitInfoView
from django.urls import path


urlpatterns = [
    path("variant_get_metrics/", VariantMetricsView.as_view(), name="variant_get_metrics"),
    path("variant_get_annotation/", VariantAnnotationView.as_view(), name="variant_get_annotation"),
    path("get_overview_stats/", OverviewDataView.as_view(), name="get_overview_stats"),
    path("overview_get_top_hits/", TopHitsView.as_view(), name="overview_get_top_hits"),
    path("trait_get_manhattan/", ManhattanView.as_view(), name="trait_get_manhattan"),
    path("trait_get_qq/", QQView.as_view(), name="trait_get_qq"),
    path("trait_get_info/", TraitInfoView.as_view(), name = "trait_get_info"),
    path("trait_get_variants/", TraitView.as_view(), name="trait_get_variants"),
    path("trait_get_chromosomeBounds/", ChromosomeBoundsView.as_view(), name="trait_get_chromosomeBounds"),
]

