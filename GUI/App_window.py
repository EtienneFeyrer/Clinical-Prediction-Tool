# -*- coding: utf-8 -*-
"""
NGS Variant Browser – Clean and optimized version with proper annotation mapping and HPO term filtering.
"""
import json
import os
import sys
import re
from dataclasses import dataclass
from collections import OrderedDict
from typing import List, Dict, Optional, Tuple, Set
import requests

from PyQt5.QtCore import QThread, pyqtSignal, Qt
from PyQt5.QtWidgets import (
    QApplication, QWidget, QHBoxLayout, QVBoxLayout, QSplitter, QTableWidget,
    QTableWidgetItem, QPushButton, QLabel, QFileDialog, QLineEdit, QFormLayout,
    QDoubleSpinBox, QMessageBox, QComboBox, QStackedWidget, QGroupBox,
    QScrollArea, QFrame, QListWidget, QListWidgetItem, QToolButton, QMenu,
    QAction, QWidgetAction, QCheckBox, QTextEdit, QSizePolicy, QTreeWidget,
    QTreeWidgetItem, QTabWidget, QTextBrowser, QDialog
)

# VCF Constants Columns
BASE_COLUMNS = [
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
    "INFO", "FORMAT", "Selector Values"
]

# Annotation Columns
ANNOTATION_COLUMNS = [
    "Gene", "Transcript ID", "Impact", "Consequence", "LOFTEE",
    "gnomAD AF", "Max Allele Freq", "SpliceAI", "GERP++",
    "PolyPhen-2", "REVEL", "CADD", "OMIM", "ML Score",
    "MANE", "cDNA Notation", "Protein Notation", "ClinSig"
]

#HPO Columns
HPO_COLUMNS = ["HPO Terms", "HPO Descriptions", "Gene-HPO Score"]

#Build columns
ALL_COLUMNS = BASE_COLUMNS + ANNOTATION_COLUMNS + HPO_COLUMNS

#define nummeric columns
NUMERIC_COLUMNS = {
    "POS", "QUAL", "gnomAD AF", "Max Allele Freq",
    "SpliceAI", "GERP++", "PolyPhen-2", "REVEL", "CADD", "ML Score", "Gene-HPO Score"
}

INTEGER_COLUMNS = {"POS"}

# HPO API Configuration (irrelevant, now local through text file)
HPO_API_BASE = "https://hpo.jax.org/api/hpo"

# Configuration for API polling
DEFAULT_API_URL = "http://127.0.0.1:5001"
POLLING_TIMEOUT = 220  # seconds
POLL_INTERVAL = 1000  # milliseconds

#filter class operators
@dataclass
class FilterRule:
    """Represents a single filter rule for table data."""
    column: str
    operator: str
    value: str = ""
    value2: str = ""
    is_numeric: bool = False
    is_automatic: bool = False

    def describe(self) -> str:
        """Return human-readable description of the filter."""
        if self.operator == "between":
            return f'{self.column} between {self.value} and {self.value2}'
        elif self.operator in ("contains", "not contains", "starts with", "ends with"):
            return f'{self.column} {self.operator} "{self.value}"'
        else:
            return f'{self.column} {self.operator} {self.value}'

# HPO filters
class HPOFilterRule(FilterRule):
    """Custom filter rule for HPO terms."""

    def __init__(self, hpo_id: str, hpo_name: str, associated_genes: Set[str]):
        super().__init__(column="Gene", operator="hpo_match", is_numeric=False)
        self.hpo_id = hpo_id
        self.hpo_name = hpo_name
        self.associated_genes = associated_genes

    def describe(self) -> str:
        return f'HPO: {self.hpo_name} ({self.hpo_id}) - {len(self.associated_genes)} genes'


class HPOManager:
    """Manager for HPO term operations using local files."""

    def __init__(self):
        self.hpo_terms = {}  # HPO_ID -> {name, definition, ...}
        self.gene_hpo_map = {}  # Gene -> set of HPO terms
        self.hpo_gene_map = {}  # HPO_ID -> set of genes
        self.loaded = False

    def load_hpo_files(self, hp_obo_path: str, phenotype_genes_path: str):
        """Load HPO data from local files."""
        print("Loading HPO data from files...")

        try:
            # Load HPO terms from .obo file
            self._parse_hp_obo(hp_obo_path)

            # Load gene-phenotype associations
            self._parse_phenotype_to_genes(phenotype_genes_path)

            self.loaded = True
            print(f"Loaded {len(self.hpo_terms)} HPO terms and {len(self.gene_hpo_map)} gene associations")

        except Exception as e:
            print(f"Error loading HPO files: {e}")
            self.loaded = False

    def _parse_hp_obo(self, file_path: str):
        """Parse the hp.obo file to extract HPO terms."""
        current_term = {}

        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()

                if line == "[Term]":
                    # Save previous term
                    if current_term and 'id' in current_term:
                        hpo_id = current_term['id']
                        self.hpo_terms[hpo_id] = {
                            'id': hpo_id,
                            'name': current_term.get('name', ''),
                            'definition': current_term.get('def', '').split('"')[1] if '"' in current_term.get('def', '') else current_term.get('def', ''),
                            'synonyms': current_term.get('synonyms', [])
                        }
                    current_term = {}

                elif line.startswith("id: HP:"):
                    current_term['id'] = line[4:]

                elif line.startswith("name: "):
                    current_term['name'] = line[6:]

                elif line.startswith("def: "):
                    current_term['def'] = line[5:]

                elif line.startswith("synonym: "):
                    if 'synonyms' not in current_term:
                        current_term['synonyms'] = []
                    # Extract synonym from quotes
                    synonym_text = line[9:]
                    if '"' in synonym_text:
                        synonym = synonym_text.split('"')[1]
                        current_term['synonyms'].append(synonym)

        # Don't forget the last term
        if current_term and 'id' in current_term:
            hpo_id = current_term['id']
            self.hpo_terms[hpo_id] = {
                'id': hpo_id,
                'name': current_term.get('name', ''),
                'definition': current_term.get('def', '').split('"')[1] if '"' in current_term.get('def', '') else current_term.get('def', ''),
                'synonyms': current_term.get('synonyms', [])
            }

    def _parse_phenotype_to_genes(self, file_path: str):
        """Parse the phenotype_to_genes.txt file."""
        with open(file_path, 'r', encoding='utf-8') as f:
            # Skip header
            header = next(f)

            for line_num, line in enumerate(f, 2):  # Start from line 2
                line = line.strip()
                if not line:
                    continue

                parts = line.split('\t')
                if len(parts) >= 4:
                    hpo_id = parts[0]
                    gene_symbol = parts[3].upper()

                    # Debug: Print first few entries
                    if line_num <= 10:
                        print(f"Line {line_num}: HPO={hpo_id}, Gene={gene_symbol}")

                    # Build gene -> HPO mapping
                    if gene_symbol not in self.gene_hpo_map:
                        self.gene_hpo_map[gene_symbol] = set()
                    self.gene_hpo_map[gene_symbol].add(hpo_id)

                    # Build HPO -> gene mapping
                    if hpo_id not in self.hpo_gene_map:
                        self.hpo_gene_map[hpo_id] = set()
                    self.hpo_gene_map[hpo_id].add(gene_symbol)
                else:
                    print(f"Skipping malformed line {line_num}: {parts}")

        print(f"Total genes loaded: {len(self.gene_hpo_map)}")
        print(f"Sample genes: {list(self.gene_hpo_map.keys())[:10]}")

    def search_hpo_terms(self, query: str, limit: int = 50) -> List[Dict]:
        """Search HPO terms by name or description."""
        if not self.loaded:
            return []

        query_lower = query.lower()
        results = []

        for hpo_id, term_data in self.hpo_terms.items():
            # Search in name
            if query_lower in term_data['name'].lower():
                results.append(term_data)
                continue

            # Search in definition
            if query_lower in term_data['definition'].lower():
                results.append(term_data)
                continue

            # Search in synonyms
            for synonym in term_data.get('synonyms', []):
                if query_lower in synonym.lower():
                    results.append(term_data)
                    break

            if len(results) >= limit:
                break

        return results[:limit]

    def get_hpo_term_details(self, hpo_id: str) -> Optional[Dict]:
        """Get detailed information for an HPO term."""
        return self.hpo_terms.get(hpo_id)

    def get_genes_for_hpo(self, hpo_id: str) -> Set[str]:
        """Get genes associated with an HPO term."""
        return self.hpo_gene_map.get(hpo_id, set())

    def get_hpo_terms_for_gene(self, gene_symbol: str) -> Set[str]:
        """Get HPO terms associated with a gene."""
        return self.gene_hpo_map.get(gene_symbol.upper(), set())

# Annotation Logic for GUI
class AnnotationWorker(QThread):
    """Worker thread for handling API communication asynchronously."""

    finished = pyqtSignal(int, int)  # successful, failed counts
    error = pyqtSignal(str)
    status_update = pyqtSignal(str)
    variant_completed = pyqtSignal(int, dict)  # row_index, annotation_data

    def __init__(self, variants_data: List[Tuple], api_url: str = DEFAULT_API_URL):
        super().__init__()
        self.variants_data = variants_data
        self.api_url = api_url
        self.successful = 0
        self.failed = 0

    def run(self):
        """Main worker execution - submit variants then poll for results."""
        try:
            self._submit_variants()
            self._poll_for_completion()
        except Exception as e:
            self.error.emit(f"Unexpected error: {str(e)}")

    def _submit_variants(self):
        """Submit all variants to the annotation service."""
        total_variants = len(self.variants_data)
        self.status_update.emit(f"Submitting {total_variants} variants...")

        self.submitted_variants = []

        for i, (row_index, chrom, pos, ref, alt) in enumerate(self.variants_data):
            variant_data = {'chrom': chrom, 'pos': int(pos), 'ref': ref, 'alt': alt}

            try:
                response = requests.post(f"{self.api_url}/submit", json=variant_data)
                variant_id = f"{chrom}:{pos}:{ref}>{alt}"

                if response.status_code == 200:
                    result = response.json()

                    if result.get('status') in ['success', 'failure']:
                        self.submitted_variants.append((row_index, variant_id))
                        status = "already annotated" if 'already annotated' in result.get('message', '') else "submitted"
                        self.status_update.emit(f"{status.title()} {i+1}/{total_variants}: {variant_id}")
                    else:
                        self._handle_submission_failure(i + 1, total_variants, result.get('message', 'Unknown error'))
                else:
                    self._handle_submission_failure(i + 1, total_variants, f"HTTP {response.status_code}")

            except Exception as e:
                self._handle_submission_failure(i + 1, total_variants, str(e))

        if not self.submitted_variants:
            self.error.emit("No variants were successfully submitted")

    def _handle_submission_failure(self, current: int, total: int, message: str):
        """Handle submission failure for a single variant."""
        self.status_update.emit(f"Failed {current}/{total}: {message}")
        self.failed += 1

    def _poll_for_completion(self):
        """Poll submitted variants until completion or timeout."""
        if not self.submitted_variants:
            return

        self.status_update.emit(f"Polling {len(self.submitted_variants)} variants...")
        completed_variants = set()

        for poll_round in range(POLLING_TIMEOUT):
            newly_completed = 0

            for row_index, variant_id in self.submitted_variants:
                if variant_id in completed_variants:
                    continue

                if self._check_variant_completion(row_index, variant_id):
                    completed_variants.add(variant_id)
                    newly_completed += 1
                    self.successful += 1

            completed_count = len(completed_variants)
            total_submitted = len(self.submitted_variants)

            self.status_update.emit(f"Poll {poll_round+1}: {completed_count}/{total_submitted} completed")

            if completed_count == total_submitted:
                break

            if poll_round > 30 and newly_completed == 0 and poll_round % 10 == 0:
                self.status_update.emit(f"Still waiting... {completed_count}/{total_submitted} completed")

            self.msleep(POLL_INTERVAL)

        self.failed += len(self.submitted_variants) - len(completed_variants)
        self.finished.emit(self.successful, self.failed)

    def _check_variant_completion(self, row_index: int, variant_id: str) -> bool:
        """Check if a specific variant is completed and emit Result if so."""
        try:
            response = requests.get(f"{self.api_url}/poll/{variant_id}")
            if response.status_code == 200:
                result = response.json()
                if result['status'] == 'completed':
                    annotation_data = result.get('annotation', {})
                    self.variant_completed.emit(row_index, annotation_data)
                    return True
        except Exception:
            pass  # Continue polling other variants
        return False

# Info Box for header inspection
class VCFHeaderInspector(QGroupBox):
    """Widget for inspecting VCF header information."""

    HEADER_PATTERN = re.compile(r'^##(INFO|FILTER|FORMAT)=<(.+)>$')

    def __init__(self, title="Header inspector (##INFO / ##FILTER / ##FORMAT)", parent=None):
        super().__init__(title, parent)
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        self._setup_ui()
        self._header_map = {"INFO": OrderedDict(), "FILTER": OrderedDict(), "FORMAT": OrderedDict()}

    def _setup_ui(self):
        """Initialize the UI components."""
        layout = QVBoxLayout(self)

        # Section and ID selection
        self.section_combo = QComboBox()
        self.section_combo.addItems(["INFO", "FILTER", "FORMAT"])

        self.id_combo = QComboBox()
        self.id_combo.setEnabled(False)

        # Text display for header details
        self.text_display = QTextEdit()
        self.text_display.setReadOnly(True)
        self.text_display.setMinimumHeight(110)
        self.text_display.setMaximumHeight(160)

        # Layout for combo boxes
        top_layout = QHBoxLayout()
        top_layout.addWidget(self.section_combo, 1)
        top_layout.addWidget(self.id_combo, 2)

        # Add to main layout
        layout.addLayout(top_layout)
        layout.addWidget(self.text_display)

        # Connect signals
        self.section_combo.currentTextChanged.connect(self._on_section_changed)
        self.id_combo.currentTextChanged.connect(self._on_id_changed)

    def set_header_lines(self, header_lines: List[str]):
        """Parse and store header lines."""
        self._header_map = {"INFO": OrderedDict(), "FILTER": OrderedDict(), "FORMAT": OrderedDict()}

        for line in header_lines:
            section, payload, raw = self._parse_header_line(line)
            if section and payload and "ID" in payload:
                self._header_map[section][payload["ID"]] = {"fields": payload, "raw": raw}

        self._populate_ids_for_section(self.section_combo.currentText())

    def _parse_header_line(self, line: str) -> Tuple[Optional[str], Optional[Dict], Optional[str]]:
        """Parse a single VCF header line."""
        match = self.HEADER_PATTERN.match(line.strip())
        if not match:
            return None, None, None

        section, inner = match.groups()
        parts = self._smart_split(inner)
        fields = OrderedDict()

        for part in parts:
            if "=" not in part:
                continue
            key, value = part.split("=", 1)
            value = value.strip()
            if len(value) >= 2 and value.startswith('"') and value.endswith('"'):
                value = value[1:-1]
            fields[key.strip()] = value

        return section, fields, line.strip()

    def _smart_split(self, text: str) -> List[str]:
        """Split text respecting quoted sections."""
        result = []
        buffer = []
        in_quotes = False

        for char in text:
            if char == '"':
                in_quotes = not in_quotes
                buffer.append(char)
            elif char == ',' and not in_quotes:
                result.append("".join(buffer).strip())
                buffer = []
            else:
                buffer.append(char)

        if buffer:
            result.append("".join(buffer).strip())
        return result

    def _on_section_changed(self, section: str):
        """Handle section combo box change."""
        self._populate_ids_for_section(section)

    def _populate_ids_for_section(self, section: str):
        """Populate ID combo box for selected section."""
        self.id_combo.blockSignals(True)
        self.id_combo.clear()

        ids = sorted(self._header_map.get(section, {}).keys())
        self.id_combo.addItems(ids)
        self.id_combo.setEnabled(bool(ids))

        self.id_combo.blockSignals(False)

        if ids:
            self._show_entry(section, ids[0])
        else:
            self.text_display.clear()

    def _on_id_changed(self, id_value: str):
        """Handle ID combo box change."""
        self._show_entry(self.section_combo.currentText(), id_value)

    def _show_entry(self, section: str, id_value: str):
        """Display detailed information for selected entry."""
        entry = self._header_map.get(section, {}).get(id_value)
        if not entry:
            self.text_display.clear()
            return

        fields = entry["fields"]
        lines = []

        # Description first
        if "Description" in fields:
            lines.append(f"Description: {fields['Description']}")

        # Standard fields
        for field in ("Number", "Type", "Source", "Version"):
            if field in fields:
                lines.append(f"{field}: {fields[field]}")

        # Other fields
        standard_fields = {"ID", "Description", "Number", "Type", "Source", "Version"}
        for key, value in fields.items():
            if key not in standard_fields:
                lines.append(f"{key}: {value}")

        lines.extend(["", f"Raw: {entry['raw']}"])
        self.text_display.setPlainText("\n".join(lines))

#HPO Terms Tab Box
class HPOFilterWidget(QGroupBox):
    """Widget for HPO term filtering."""

    def __init__(self, hpo_manager: HPOManager, filter_manager, parent=None):
        super().__init__("HPO Term Filters", parent)
        self.hpo_manager = hpo_manager
        self.filter_manager = filter_manager
        self.selected_hpo_terms = {}  # HPO_ID -> HPO data
        self._setup_ui()

    def _setup_ui(self):
        layout = QVBoxLayout(self)

        # HPO status label - this should be in the HPO tab
        self.hpo_status_label = QLabel("HPO files not loaded. Click 'Load HPO Files' to enable filtering.")
        self.hpo_status_label.setStyleSheet("color: orange; font-weight: bold;")
        layout.addWidget(self.hpo_status_label)

        # Add HPO matching button
        self.match_hpo_btn = QPushButton("Match HPO Terms to Genes")
        self.match_hpo_btn.setEnabled(False)
        self.match_hpo_btn.setStyleSheet("font-weight: bold; padding: 8px;")
        layout.addWidget(self.match_hpo_btn)

        # Search section
        search_group = QGroupBox("Search HPO Terms")
        search_layout = QVBoxLayout(search_group)

        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText("Search HPO terms (e.g., 'seizure', 'intellectual disability')...")
        self.search_button = QPushButton("Search")

        # Initially disable until HPO files are loaded
        self.search_input.setEnabled(False)
        self.search_button.setEnabled(False)

        search_top = QHBoxLayout()
        search_top.addWidget(self.search_input)
        search_top.addWidget(self.search_button)
        search_layout.addLayout(search_top)

        # Results tree
        self.results_tree = QTreeWidget()
        self.results_tree.setHeaderLabels(["HPO ID", "Name", "Description"])
        self.results_tree.setRootIsDecorated(False)
        self.results_tree.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        search_layout.addWidget(self.results_tree)

        # Action buttons
        button_layout = QHBoxLayout()
        self.add_filter_btn = QPushButton("Add Selected as Filter")
        self.view_genes_btn = QPushButton("View Associated Genes")
        button_layout.addWidget(self.add_filter_btn)
        button_layout.addWidget(self.view_genes_btn)
        search_layout.addLayout(button_layout)


        layout.addWidget(search_group)

        # Active HPO filters
        active_group = QGroupBox("Active HPO Filters")
        active_layout = QVBoxLayout(active_group)

        self.active_tree = QTreeWidget()
        self.active_tree.setHeaderLabels(["HPO ID", "Name", "Matching Variants"])
        self.active_tree.setRootIsDecorated(False)

        self.clear_hpo_btn = QPushButton("Clear All HPO Filters")

        active_layout.addWidget(self.active_tree)
        active_layout.addWidget(self.clear_hpo_btn)
        layout.addWidget(active_group)

        # Connect signals
        self.search_button.clicked.connect(self._search_hpo_terms)
        self.search_input.returnPressed.connect(self._search_hpo_terms)
        self.add_filter_btn.clicked.connect(self._add_hpo_filter)
        self.view_genes_btn.clicked.connect(self._view_associated_genes)
        self.clear_hpo_btn.clicked.connect(self._clear_hpo_filters)
        self.active_tree.itemDoubleClicked.connect(self._remove_hpo_filter)
        self.match_hpo_btn.clicked.connect(self._match_hpo_terms_to_genes)

    def _search_hpo_terms(self):
        query = self.search_input.text().strip()
        if not query:
            return

        self.results_tree.clear()
        self.search_button.setEnabled(False)
        self.search_button.setText("Searching...")

        try:
            terms = self.hpo_manager.search_hpo_terms(query)

            for term in terms:
                item = QTreeWidgetItem([
                    term.get('id', ''),
                    term.get('name', ''),
                    term.get('definition', '')[:100] + '...' if len(term.get('definition', '')) > 100
                    else term.get('definition', '')
                ])
                item.setData(0, Qt.UserRole, term)
                self.results_tree.addTopLevelItem(item)

        except Exception as e:
            QMessageBox.warning(self, "Search Error", f"Error searching HPO terms: {e}")
        finally:
            self.search_button.setEnabled(True)
            self.search_button.setText("Search")

    def _add_hpo_filter(self):
        selected_items = self.results_tree.selectedItems()
        if not selected_items:
            QMessageBox.information(self, "No Selection", "Please select an HPO term to add as filter.")
            return

        for item in selected_items:
            term_data = item.data(0, Qt.UserRole)
            if term_data:
                hpo_id = term_data.get('id', '')
                if hpo_id not in self.selected_hpo_terms:
                    self.selected_hpo_terms[hpo_id] = term_data
                    self._apply_hpo_filter(hpo_id, term_data)

        self._update_active_filters_display()

    def _apply_hpo_filter(self, hpo_id: str, term_data: Dict):
        """Apply HPO filter by finding genes associated with the HPO term."""
        # Get genes for this HPO term
        associated_genes = self.hpo_manager.get_genes_for_hpo(hpo_id)

        print(f"Applying HPO filter for {hpo_id}: {term_data.get('name', '')}")
        print(f"Associated genes ({len(associated_genes)}): {list(associated_genes)[:10]}...")

        if not associated_genes:
            QMessageBox.information(
                self, "No Genes Found",
                f"No genes found associated with {term_data.get('name', hpo_id)}"
            )
            return

        # Create a custom filter rule for HPO
        rule = HPOFilterRule(
            hpo_id=hpo_id,
            hpo_name=term_data.get('name', ''),
            associated_genes=associated_genes
        )

        self.filter_manager.add_filter(rule)
# Genes associated
    def _view_associated_genes(self):
        selected_items = self.results_tree.selectedItems()
        if not selected_items:
            QMessageBox.information(self, "No Selection", "Please select an HPO term to view genes.")
            return

        item = selected_items[0]
        term_data = item.data(0, Qt.UserRole)
        if not term_data:
            return

        hpo_id = term_data.get('id', '')
        genes = self.hpo_manager.get_genes_for_hpo(hpo_id)

        # Show genes in a dialog
        dialog = QDialog(self)
        dialog.setWindowTitle(f"Genes for {term_data.get('name', hpo_id)}")
        dialog.setMinimumSize(400, 300)

        layout = QVBoxLayout(dialog)

        text_browser = QTextBrowser()
        if genes:
            gene_list = sorted(genes)
            text_browser.setHtml(f"""
            <h3>{term_data.get('name', hpo_id)} ({hpo_id})</h3>
            <p><strong>Description:</strong> {term_data.get('definition', 'No description available')}</p>
            <p><strong>Associated Genes ({len(genes)}):</strong></p>
            <p>{', '.join(gene_list)}</p>
            """)
        else:
            text_browser.setHtml("<p>No genes found for this HPO term.</p>")

        layout.addWidget(text_browser)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        layout.addWidget(close_btn)

        dialog.exec_()

    def _clear_hpo_filters(self):
        # Remove all HPO filters
        hpo_rules_to_remove = [rule for rule in self.filter_manager.filters
                               if isinstance(rule, HPOFilterRule)]

        for rule in hpo_rules_to_remove:
            self.filter_manager.remove_filter(rule)

        self.selected_hpo_terms.clear()
        self._update_active_filters_display()

    def _remove_hpo_filter(self, item):
        hpo_id = item.text(0)
        if hpo_id in self.selected_hpo_terms:
            # Find and remove the corresponding filter rule
            rule_to_remove = None
            for rule in self.filter_manager.filters:
                if isinstance(rule, HPOFilterRule) and rule.hpo_id == hpo_id:
                    rule_to_remove = rule
                    break

            if rule_to_remove:
                self.filter_manager.remove_filter(rule_to_remove)
                del self.selected_hpo_terms[hpo_id]
                self._update_active_filters_display()

    def _update_active_filters_display(self):
        self.active_tree.clear()

        for hpo_id, term_data in self.selected_hpo_terms.items():
            # Count matching variants
            matching_count = self._count_matching_variants(hpo_id)

            item = QTreeWidgetItem([
                hpo_id,
                term_data.get('name', ''),
                str(matching_count)
            ])
            self.active_tree.addTopLevelItem(item)

    def _count_matching_variants(self, hpo_id: str) -> int:
        """Count variants that match the HPO filter."""
        if not hasattr(self.filter_manager, 'table'):
            return 0

        count = 0
        gene_col_idx = ALL_COLUMNS.index("Gene")
        associated_genes = self.hpo_manager.get_genes_for_hpo(hpo_id)

        for row in range(self.filter_manager.table.rowCount()):
            if self.filter_manager.table.isRowHidden(row):
                continue

            gene_item = self.filter_manager.table.item(row, gene_col_idx)
            gene_value = gene_item.text().upper() if gene_item else ""

            if gene_value in associated_genes:
                count += 1

        return count
    def set_hpo_loaded(self, loaded: bool):
        self.search_input.setEnabled(loaded)
        self.search_button.setEnabled(loaded)
        self.add_filter_btn.setEnabled(loaded)
        self.view_genes_btn.setEnabled(loaded)
        self.match_hpo_btn.setEnabled(loaded)  # Add this line

        if loaded:
            terms_count = len(self.hpo_manager.hpo_terms) if hasattr(self.hpo_manager, 'hpo_terms') else 0
            genes_count = len(self.hpo_manager.gene_hpo_map) if hasattr(self.hpo_manager, 'gene_hpo_map') else 0
            self.hpo_status_label.setText(f"HPO data loaded: {terms_count} terms, {genes_count} genes")
            self.hpo_status_label.setStyleSheet("color: green; font-weight: bold;")
        else:
            self.hpo_status_label.setText("HPO files not loaded. Click 'Load HPO Files' to enable filtering.")
            self.hpo_status_label.setStyleSheet("color: orange; font-weight: bold;")


    def _match_hpo_terms_to_genes(self):
        """Match HPO terms to genes in the current table."""
        if not self.hpo_manager.loaded:
            QMessageBox.warning(self, "HPO Not Loaded", "Please load HPO files first.")
            return

        if not hasattr(self.filter_manager, 'table') or self.filter_manager.table.rowCount() == 0:
            QMessageBox.information(self, "No Data", "No variants loaded in table.")
            return

        # Show progress dialog
        progress_dialog = QMessageBox(self)
        progress_dialog.setWindowTitle("Matching HPO Terms")
        progress_dialog.setText("Matching HPO terms to genes... Please wait.")
        progress_dialog.setStandardButtons(QMessageBox.NoButton)
        progress_dialog.setModal(True)
        progress_dialog.show()
        QApplication.processEvents()

        try:
            # Direct access to the table and HPO population
            matched_count = self._populate_hpo_data_directly()

            progress_dialog.hide()
            progress_dialog.deleteLater()

            if matched_count > 0:
                QMessageBox.information(
                    self, "HPO Matching Complete",
                    f"Successfully matched HPO terms for {matched_count} variants with gene information."
                )
            else:
                QMessageBox.information(
                    self, "No Matches Found",
                    "No variants with gene information found. Please:\n"
                    "1. Load a VCF file with gene annotations, or\n"
                    "2. Use 'Annotate via API' to get gene information, or\n"
                    "3. Manually add test genes for demonstration."
                )

        except Exception as e:
            progress_dialog.hide()
            progress_dialog.deleteLater()
            QMessageBox.critical(self, "Error", f"Error during HPO matching: {str(e)}")

    def _populate_hpo_data_directly(self):
        """Populate HPO data directly from the filter widget."""
        if not self.hpo_manager.loaded:
            return 0

        table = self.filter_manager.table
        gene_col_idx = ALL_COLUMNS.index("Gene")
        hpo_terms_col = ALL_COLUMNS.index("HPO Terms")
        hpo_desc_col = ALL_COLUMNS.index("HPO Descriptions")
        hpo_score_col = ALL_COLUMNS.index("Gene-HPO Score")

        populated_count = 0

        for row in range(table.rowCount()):
            gene_item = table.item(row, gene_col_idx)
            if gene_item and gene_item.text().strip():
                gene_symbol = gene_item.text().strip()

                # Get HPO terms for this gene
                hpo_terms = self.hpo_manager.get_hpo_terms_for_gene(gene_symbol)

                if hpo_terms:
                    # Get details for HPO terms
                    term_names = []
                    term_descriptions = []

                    for hpo_id in list(hpo_terms)[:5]:  # Limit to first 5 terms
                        term_data = self.hpo_manager.get_hpo_term_details(hpo_id)
                        if term_data:
                            term_names.append(f"{hpo_id}:{term_data.get('name', '')}")
                            term_descriptions.append(term_data.get('definition', '')[:50] + '...' if len(term_data.get('definition', '')) > 50 else term_data.get('definition', ''))

                    # Calculate a simple score based on number of associated HPO terms
                    hpo_score = min(len(hpo_terms) / 10.0, 1.0)  # Normalize to 0-1 range

                    # Update table cells
                    table.setItem(row, hpo_terms_col, QTableWidgetItem("; ".join(term_names)))
                    table.setItem(row, hpo_desc_col, QTableWidgetItem("; ".join(term_descriptions)))
                    table.setItem(row, hpo_score_col, QTableWidgetItem(f"{hpo_score:.3f}"))
                    populated_count += 1
                else:
                    # Set empty values
                    table.setItem(row, hpo_terms_col, QTableWidgetItem(""))
                    table.setItem(row, hpo_desc_col, QTableWidgetItem(""))
                    table.setItem(row, hpo_score_col, QTableWidgetItem(""))

        print(f"Populated HPO data for {populated_count} rows")
        return populated_count

# Handler for the incoming Annoations
class AnnotationHandler:
    """Handles annotation data extraction and formatting."""
    @staticmethod
    def extract_annotation_data(annotation_data: dict) -> Dict[str, str]:
            """Extract all annotation data for table display."""
            return {
                    'Gene': AnnotationHandler._get_direct_value(annotation_data, ['gene']),
                    'Transcript ID': AnnotationHandler._get_transcript_value(annotation_data, 'transcript_id'),
                    'Impact': AnnotationHandler._get_transcript_value(annotation_data, 'impact'),
                    'Consequence': AnnotationHandler._get_direct_value(annotation_data, ['Most Severe Consequence']),
                    'LOFTEE': AnnotationHandler._get_transcript_value(annotation_data, 'LOFTEE'),
                    'gnomAD AF': AnnotationHandler._get_direct_value(annotation_data, ['gnomAD AF']),
                    'Max Allele Freq': AnnotationHandler._get_direct_value(annotation_data, ['max_allele_freq']),
                    'SpliceAI': AnnotationHandler._get_transcript_value(annotation_data, 'Splice_AI'),
                    'GERP++': AnnotationHandler._get_transcript_value(annotation_data, 'GERP'),
                    'PolyPhen-2': AnnotationHandler._get_transcript_value(annotation_data, 'polyphen'),
                    'REVEL': AnnotationHandler._get_transcript_value(annotation_data, 'REVEL'),
                    'CADD': AnnotationHandler._get_direct_value(annotation_data, ['CADD']),
                    'OMIM': AnnotationHandler._get_direct_value(annotation_data, ['OMIM']),
                    'ML Score': AnnotationHandler._get_direct_value(annotation_data, ['pathogenicity_score']),
                    'MANE': AnnotationHandler._get_transcript_value(annotation_data, 'Mane'),
                    'cDNA Notation': AnnotationHandler._get_transcript_value(annotation_data, 'cDNA_notation'),
                    'Protein Notation': AnnotationHandler._get_transcript_value(annotation_data, 'protein_notation'),
                    'ClinSig': AnnotationHandler._get_direct_value(annotation_data, ['CLINSIG'])  # Add this line
                }

    @staticmethod
    def _get_direct_value(data: dict, keys: List[str]) -> str:
        """Get value directly from annotation data."""
        for key in keys:
            if key in data and data[key] is not None:
                return str(data[key]).strip()
            # Case-insensitive fallback
            for data_key in data.keys():
                if data_key.lower() == key.lower() and data[data_key] is not None:
                    return str(data[data_key]).strip()
        return ""

    @staticmethod
    def _get_transcript_value(data: dict, field: str) -> str:
        """Get value from first transcript."""
        transcripts = data.get('transcript_consequences', [])
        if transcripts and field in transcripts[0]:
            value = transcripts[0][field]
            if value is not None:
                # Special handling for boolean MANE field
                if field == 'Mane':
                    return "Yes" if value else "No"
                return str(value).strip()
        return ""


    @staticmethod
    def format_value(column_name: str, value: str) -> str:
        """Format value for display based on column type."""
        if not value:
            return ""

        value_str = str(value).strip()

        # ClinSig formatting
        if column_name == 'ClinSig':
            if isinstance(value_str, str):
                # Remove brackets if present and clean up
                cleaned = value_str.strip("[]'\"").replace("_", " ")
                return cleaned.title() if cleaned else ""
            return value_str

        # Numeric scores
        if column_name in ['REVEL', 'SpliceAI', 'GERP++', 'CADD', 'ML Score', 'PolyPhen-2']:
            try:
                num_val = float(value_str)
                return f"{num_val:.3f}" if num_val != int(num_val) else str(int(num_val))
            except ValueError:
                pass

        # Frequencies
        elif column_name in ['gnomAD AF', 'Max Allele Freq']:
            try:
                freq_val = float(value_str)
                if 0 < freq_val < 0.01:
                    return f"{freq_val*100:.4f}%"
                elif freq_val < 1:
                    return f"{freq_val:.6f}"
            except ValueError:
                pass

        # Truncate long values
        if len(value_str) > 100:
            return value_str[:97] + "..."

        return value_str

# In Table filtering
class FilterManager:
    """Manages table filtering operations."""

    def __init__(self, table: QTableWidget):
        self.table = table
        self.filters: List[FilterRule] = []

    def add_filter(self, rule: FilterRule):
        """Add a new filter rule."""
        self.filters.append(rule)
        self.apply_filters()

    def remove_filter(self, rule: FilterRule):
        """Remove a filter rule."""
        try:
            self.filters.remove(rule)
            self.apply_filters()
        except ValueError:
            pass

    def clear_filters(self):
        """Remove all filters."""
        self.filters.clear()
        for row in range(self.table.rowCount()):
            self.table.setRowHidden(row, False)

    def clear_automatic_filters(self):  #
        """Remove only automatic filters."""
        self.filters = [rule for rule in self.filters if not rule.is_automatic]
        self.apply_filters()

    def has_automatic_filters(self) -> bool:
        """Check if there are any automatic filters active."""
        return any(rule.is_automatic for rule in self.filters)

    def apply_filters(self):
        """Apply all active filters to the table."""
        if not self.filters:
            self.clear_filters()
            return

        for row in range(self.table.rowCount()):
            row_visible = all(self._match_rule(rule, row) for rule in self.filters)
            self.table.setRowHidden(row, not row_visible)

    def _match_rule(self, rule: FilterRule, row: int) -> bool:
        """Check if a row matches a specific filter rule."""
        # Handle HPO filter rules
        if isinstance(rule, HPOFilterRule):
            return self._match_hpo_rule(rule, row)

        col_idx = ALL_COLUMNS.index(rule.column)
        item = self.table.item(row, col_idx)
        text_value = item.text() if item else ""

        if rule.is_numeric:
            return self._match_numeric_rule(rule, text_value)
        else:
            return self._match_text_rule(rule, text_value)

    def _match_hpo_rule(self, rule: HPOFilterRule, row: int) -> bool:
        """Check if a row matches an HPO filter rule."""
        # Get the gene from the row
        gene_col_idx = ALL_COLUMNS.index("Gene")
        gene_item = self.table.item(row, gene_col_idx)
        gene_value = gene_item.text().strip().upper() if gene_item else ""

        if not gene_value:
            return False

        # Debug: Print gene matching info
        print(f"Checking gene '{gene_value}' against HPO term {rule.hpo_id}")
        print(f"Associated genes: {list(rule.associated_genes)[:10]}...")  # Show first 10 genes

        # Check if this gene is associated with the HPO term
        match_found = gene_value in rule.associated_genes
        if match_found:
            print(f"MATCH FOUND: {gene_value}")

        return match_found

    def _match_numeric_rule(self, rule: FilterRule, text_value: str) -> bool:
        """Match numeric filter rule."""
        try:
            value = float(text_value.strip()) if text_value.strip() else None
        except ValueError:
            return False

        if value is None:
            return False

        if rule.operator == "between":
            return float(rule.value) <= value <= float(rule.value2)

        target = float(rule.value)
        operators = {
            ">": lambda v, t: v > t,
            ">=": lambda v, t: v >= t,
            "<": lambda v, t: v < t,
            "<=": lambda v, t: v <= t,
            "=": lambda v, t: v == t,
            "!=": lambda v, t: v != t
        }
        return operators.get(rule.operator, lambda v, t: True)(value, target)

    def _match_text_rule(self, rule: FilterRule, text_value: str) -> bool:
        """Match text filter rule."""
        text = text_value.lower()
        target = rule.value.lower()

        operators = {
            "contains": lambda t, tgt: tgt in t,
            "not contains": lambda t, tgt: tgt not in t,
            "starts with": lambda t, tgt: t.startswith(tgt),
            "ends with": lambda t, tgt: t.endswith(tgt),
            "=": lambda t, tgt: t == tgt,
            "!=": lambda t, tgt: t != tgt
        }
        return operators.get(rule.operator, lambda t, tgt: True)(text, target)

# The main Window
class AppWindow(QWidget):
    """Main application window."""

    def __init__(self):
        super().__init__()
        self.hpo_manager = HPOManager()
        self._setup_window()
        self._init_data()
        self._setup_ui()
        self._connect_signals()
        self._initialize_components()

    def _setup_window(self):
        """Configure main window properties."""
        self.setWindowTitle("NGS Variant Browser")
        self.resize(1450, 900)

    def _init_data(self):
        """Initialize data structures."""
        self.headers: List[str] = []
        self.variants: List[List[str]] = []
        self.dynamic_format_columns: List[str] = []

        # Set default column visibility Change on requirement
        self.column_visibility: Dict[str, bool] = {}
        for col in ALL_COLUMNS:
            if col in ["FORMAT", "FILTER", "INFO", "Selector Values"]:
                self.column_visibility[col] = False  # Hide these columns by default!!!
            else:
                self.column_visibility[col] = True   # Show all other columns by default

        # Components
        self.annotation_worker = None
        self.filter_manager = None
        self.api_url = DEFAULT_API_URL

    def _setup_ui(self):
        """Create and arrange UI components."""
        # Main table
        self.table = QTableWidget()
        self.table.setColumnCount(len(ALL_COLUMNS))
        self.table.setHorizontalHeaderLabels(ALL_COLUMNS)
        self.table.setSortingEnabled(True)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)

        # Control buttons
        self.load_button = QPushButton("Load VCF File")
        self.save_json_button = QPushButton("Save JSON File")
        self.annotate_button = QPushButton("Annotate via API")
        self.status_label = QLabel("No file loaded")
        self.load_hpo_button = QPushButton("Load HPO Files")

        self.auto_filter_button = QPushButton("Auto-Filter: ON")
        self.auto_filter_button.setCheckable(True)
        self.auto_filter_button.setChecked(True)
        self.auto_filter_button.setStyleSheet("QPushButton:checked { background-color: #4CAF50; }")


    # Create tabbed interface for filters
        self.filter_tabs = QTabWidget()

        # Regular filters tab
        self.regular_filter_widget = QWidget()
        self._setup_regular_filters(self.regular_filter_widget)
        self.filter_tabs.addTab(self.regular_filter_widget, "Regular Filters")

        # Header inspector
        self.header_inspector = VCFHeaderInspector(parent=self)

        # Layout
        self._arrange_layout()

    def _setup_regular_filters(self, parent_widget):
        """Setup the regular filter UI in the provided parent widget."""
        layout = QVBoxLayout(parent_widget)

        # Filter components
        # Column and operator selection
        self.column_combo = QComboBox()
        self.column_combo.addItems(ALL_COLUMNS)
        self.operator_combo = QComboBox()

        # Input stack for different value types
        self.input_stack = QStackedWidget()

        # Text input
        self.input_text = QLineEdit()
        self.input_text.setPlaceholderText("Enter text…")
        text_widget = QWidget()
        text_layout = QVBoxLayout(text_widget)
        text_layout.setContentsMargins(0, 0, 0, 0)
        text_layout.addWidget(self.input_text)
        self.input_stack.addWidget(text_widget)

        # Numeric input
        self.input_num = QDoubleSpinBox()
        self.input_num.setRange(-1e12, 1e12)
        self.input_num.setDecimals(6)
        num_widget = QWidget()
        num_layout = QVBoxLayout(num_widget)
        num_layout.setContentsMargins(0, 0, 0, 0)
        num_layout.addWidget(self.input_num)
        self.input_stack.addWidget(num_widget)

        # Range input
        self.input_num_min = QDoubleSpinBox()
        self.input_num_min.setRange(-1e12, 1e12)
        self.input_num_min.setDecimals(6)
        self.input_num_max = QDoubleSpinBox()
        self.input_num_max.setRange(-1e12, 1e12)
        self.input_num_max.setDecimals(6)
        range_widget = QWidget()
        range_layout = QHBoxLayout(range_widget)
        range_layout.setContentsMargins(0, 0, 0, 0)
        range_layout.addWidget(self.input_num_min)
        range_layout.addWidget(QLabel("and"))
        range_layout.addWidget(self.input_num_max)
        self.input_stack.addWidget(range_widget)

        # Filter control buttons
        self.add_filter_btn = QPushButton("Add filter")
        self.reset_filters_btn = QPushButton("Reset all filters")

        # Values list for quick filtering
        self.values_search = QLineEdit()
        self.values_search.setPlaceholderText("Filter values")
        self.values_list = QListWidget()
        self.values_list.setSelectionMode(QListWidget.SingleSelection)
        self.add_value_btn = QPushButton('Add "=" filter for selected value')

        # Active filters display
        self.filters_area = QScrollArea()
        self.filters_area.setWidgetResizable(True)
        self.filters_host = QWidget()
        self.filters_vbox = QVBoxLayout(self.filters_host)
        self.filters_vbox.setContentsMargins(8, 8, 8, 8)
        self.filters_vbox.addStretch(1)
        self.filters_area.setWidget(self.filters_host)

        # Column visibility
        self.columns_button = QToolButton()
        self.columns_button.setText("Choose columns")
        self.columns_button.setPopupMode(QToolButton.InstantPopup)
        self.columns_menu = QMenu(self)
        self.columns_button.setMenu(self.columns_menu)

        # Add filter group
        add_group = QGroupBox("Add Filter")
        add_layout = QVBoxLayout(add_group)
        form_layout = QFormLayout()
        form_layout.addRow("Column", self.column_combo)
        form_layout.addRow("Operator", self.operator_combo)
        form_layout.addRow("Value", self.input_stack)
        add_layout.addLayout(form_layout)
        add_layout.addWidget(self.add_filter_btn)

        # Values group
        values_group = QGroupBox("Distinct values in column (from visible rows)")
        values_layout = QVBoxLayout(values_group)
        values_layout.addWidget(self.values_search)
        values_layout.addWidget(self.values_list, 1)
        values_layout.addWidget(self.add_value_btn)

        # Active filters group
        active_group = QGroupBox("Active Filters (AND)")
        active_layout = QVBoxLayout(active_group)
        active_layout.addWidget(self.filters_area)

        filter_buttons_layout = QHBoxLayout()
        filter_buttons_layout.addWidget(self.reset_filters_btn)
        filter_buttons_layout.addWidget(self.auto_filter_button)  # Add this line
        active_layout.addLayout(filter_buttons_layout)  # Add this line instead of adding reset_filters_btn directly


    # Column visibility group
        cols_group = QGroupBox("Columns shown")
        cols_layout = QVBoxLayout(cols_group)
        cols_layout.addWidget(self.columns_button)

        # Assemble regular filter layout
        layout.addWidget(add_group)
        layout.addWidget(values_group, 1)
        layout.addWidget(active_group)
        layout.addWidget(cols_group)

    def _arrange_layout(self):
        """Arrange all UI components in the main layout."""
        # Top bar
        top_bar = QHBoxLayout()
        top_bar.addWidget(self.load_button)
        top_bar.addWidget(self.save_json_button)
        top_bar.addWidget(self.annotate_button)
        top_bar.addWidget(self.load_hpo_button)
        top_bar.addStretch(1)
        top_bar.addWidget(self.status_label)

        # Left panel (table)
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.setContentsMargins(8, 8, 8, 8)
        left_layout.addLayout(top_bar)
        left_layout.addWidget(self.table, 1)

        # Right panel (filters and controls)
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        right_layout.setContentsMargins(8, 8, 8, 8)

        # Add the tabbed filter interface
        right_layout.addWidget(self.filter_tabs, 1)
        right_layout.addWidget(self.header_inspector)
        right_layout.addStretch(1)

        # Main splitter
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(left_panel)
        splitter.addWidget(right_panel)
        splitter.setStretchFactor(0, 5)
        splitter.setStretchFactor(1, 2)

        main_layout = QHBoxLayout(self)
        main_layout.addWidget(splitter)

    def _connect_signals(self):
        """Connect all signal handlers."""
        # Buttons
        self.load_button.clicked.connect(self._load_vcf)
        self.save_json_button.clicked.connect(self._save_json)
        self.annotate_button.clicked.connect(self._annotate_variants)
        self.load_hpo_button.clicked.connect(self._load_hpo_files)

        # Filter controls
        self.column_combo.currentTextChanged.connect(self._on_column_changed)
        self.operator_combo.currentTextChanged.connect(self._update_input_widget)
        self.add_filter_btn.clicked.connect(self._add_filter)
        self.reset_filters_btn.clicked.connect(self._reset_filters)
        self.auto_filter_button.clicked.connect(self._toggle_auto_filter)

        # Values list
        self.values_search.textChanged.connect(self._filter_values_list)
        self.add_value_btn.clicked.connect(self._add_value_filter)
        self.values_list.itemDoubleClicked.connect(self._add_value_filter)

    def _initialize_components(self):
        """Initialize component states."""
        self.filter_manager = FilterManager(self.table)

        # HPO filters tab - create after filter_manager is initialized
        self.hpo_filter_widget = HPOFilterWidget(self.hpo_manager, self.filter_manager)
        self.filter_tabs.addTab(self.hpo_filter_widget, "HPO Filters")

        self._build_column_menu()
        self._update_operator_options()
        self._update_values_list()
        self._apply_column_visibility()

    # Event Handlers
    def _load_vcf(self):
        """Load VCF file dialog and processing."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open VCF", "", "VCF files (*.vcf *.vcf.gz *.txt);;All files (*)"
        )

        if not file_path:
            return

        try:
            self.headers, self.variants = VCFReader.read_vcf(file_path)
            self._populate_table()
            self._debug_gene_column()
            self._update_header_inspector(file_path)
            self.filter_manager.clear_filters()
            self.status_label.setText(f"Loaded {len(self.variants)} variants")

        except Exception as e:
            QMessageBox.critical(self, "Load Error", f"Could not load file:\n{type(e).__name__}: {e}")

    def _save_json(self):
        """Save visible table data to JSON file."""
        if self.table.rowCount() == 0:
            QMessageBox.information(self, "Nothing to save", "No variants to export.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save JSON", "variants.json", "JSON (*.json)"
        )

        if not file_path:
            return

        try:
            data = []
            for row in range(self.table.rowCount()):
                if self.table.isRowHidden(row):
                    continue

                record = {}
                for col, column_name in enumerate(ALL_COLUMNS):
                    item = self.table.item(row, col)
                    record[column_name] = item.text() if item else ""
                data.append(record)

            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2, ensure_ascii=False)

            self.status_label.setText(f"Saved {len(data)} visible rows to {os.path.basename(file_path)}")

        except Exception as e:
            QMessageBox.critical(self, "Save Error", f"Could not save JSON:\n{type(e).__name__}: {e}")

    def _annotate_variants(self):
        """Start annotation process for visible variants."""
        if not self.variants:
            QMessageBox.information(self, "No data", "Please load a VCF file first.")
            return

        variants_to_annotate = self._collect_variants_for_annotation()

        if not variants_to_annotate:
            QMessageBox.information(self, "No variants", "No valid variants found to annotate.")
            return

        # Confirm with user
        reply = QMessageBox.question(
            self, "Confirm Annotation",
            f"Found {len(variants_to_annotate)} variant(s) to annotate. Continue?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes
        )

        if reply != QMessageBox.Yes:
            return

        self._start_annotation_worker(variants_to_annotate)

    def _collect_variants_for_annotation(self) -> List[Tuple]:
        """Collect visible variants for annotation."""
        variants = []
        col_indices = {col: ALL_COLUMNS.index(col) for col in ["#CHROM", "POS", "REF", "ALT"]}

        for row in range(self.table.rowCount()):
            if self.table.isRowHidden(row):
                continue

            # Get variant data
            variant_data = {}
            for col_name, col_idx in col_indices.items():
                item = self.table.item(row, col_idx)
                if not item or not item.text().strip():
                    break
                variant_data[col_name] = item.text().strip()
            else:
                # Handle multiple ALT alleles
                for alt_allele in variant_data["ALT"].split(','):
                    alt_allele = alt_allele.strip()
                    if alt_allele:
                        variants.append((
                            row, variant_data["#CHROM"], variant_data["POS"],
                            variant_data["REF"], alt_allele
                        ))

        return variants

    def _start_annotation_worker(self, variants_to_annotate: List[Tuple]):
        """Initialize and start annotation worker thread."""
        self.annotate_button.setEnabled(False)
        self.annotate_button.setText("Annotating...")
        self.status_label.setText("Starting annotation...")

        try:
            self.annotation_worker = AnnotationWorker(variants_to_annotate, self.api_url)
            self.annotation_worker.finished.connect(self._on_annotation_finished)
            self.annotation_worker.error.connect(self._on_annotation_error)
            self.annotation_worker.status_update.connect(self._on_annotation_status)
            self.annotation_worker.variant_completed.connect(self._on_variant_completed)
            self.annotation_worker.start()

        except Exception as e:
            self._reset_annotation_ui()
            QMessageBox.critical(self, "Annotation Error", f"Failed to start annotation:\n{e}")

    def _reset_annotation_ui(self):
        """Reset annotation UI to default state."""
        self.annotate_button.setEnabled(True)
        self.annotate_button.setText("Annotate via API")

    # Annotation worker event handlers
    def _on_annotation_finished(self, successful: int, failed: int):
        """Handle annotation completion."""
        self._reset_annotation_ui()
        total = successful + failed

        if successful > 0:
            self.status_label.setText(f"Annotation completed: {successful}/{total} successful")
            QMessageBox.information(
                self, "Annotation Complete",
                f"Annotation completed!\n\nSuccessful: {successful}\nFailed: {failed}"
            )
        else:
            self.status_label.setText("Annotation failed - no variants completed")
            QMessageBox.warning(
                self, "Annotation Failed",
                f"No variants were successfully annotated.\nFailed: {failed}"
            )

        self.annotation_worker = None

    def _on_annotation_error(self, message: str):
        """Handle annotation error."""
        self._reset_annotation_ui()
        self.status_label.setText("Annotation failed")
        QMessageBox.critical(self, "Annotation Error", message)
        self.annotation_worker = None

    def _on_annotation_status(self, text: str):
        """Handle annotation status updates."""
        self.status_label.setText(text)

    def _on_variant_completed(self, row_index: int, annotation_data: dict):
        """Handle individual variant completion."""
        try:
            annotation_mapping = AnnotationHandler.extract_annotation_data(annotation_data)

            for column_name, value in annotation_mapping.items():
                if column_name in ANNOTATION_COLUMNS and value:
                    col_idx = ALL_COLUMNS.index(column_name)
                    formatted_value = AnnotationHandler.format_value(column_name, value)
                    self.table.setItem(row_index, col_idx, QTableWidgetItem(formatted_value))

            # Update HPO data for the gene
            gene_value = annotation_mapping.get('Gene', '')
            if gene_value:
                self._populate_hpo_data(row_index, gene_value)

        except Exception as e:
            print(f"Error updating row {row_index}: {e}")

    # Filter management
    def _on_column_changed(self):
        """Handle column selection change in filter UI."""
        self._update_operator_options()
        self._update_values_list()

    def _update_operator_options(self):
        """Update operator options based on selected column."""
        column = self.column_combo.currentText()
        self.operator_combo.blockSignals(True)
        self.operator_combo.clear()

        if column in NUMERIC_COLUMNS:
            self.operator_combo.addItems(["=", "!=", ">", ">=", "<", "<=", "between"])
            self.operator_combo.setCurrentText("between")
        else:
            self.operator_combo.addItems(["contains", "not contains", "=", "!=", "starts with", "ends with"])
            self.operator_combo.setCurrentText("contains")

        self.operator_combo.blockSignals(False)
        self._update_input_widget()

    def _update_input_widget(self):
        """Update input widget based on selected operator."""
        operator = self.operator_combo.currentText()
        column = self.column_combo.currentText()

        if operator == "between":
            self.input_stack.setCurrentIndex(2)
        elif operator in (">", ">=", "<", "<=", "=", "!=") and column in NUMERIC_COLUMNS:
            self.input_stack.setCurrentIndex(1)
        else:
            self.input_stack.setCurrentIndex(0)

    def _add_filter(self):
        """Add filter from UI inputs."""
        column = self.column_combo.currentText()
        operator = self.operator_combo.currentText()
        is_numeric = column in NUMERIC_COLUMNS

        try:
            if operator == "between":
                rule = FilterRule(
                    column=column, operator=operator,
                    value=str(self.input_num_min.value()),
                    value2=str(self.input_num_max.value()),
                    is_numeric=True
                )
            elif is_numeric and operator in (">", ">=", "<", "<=", "=", "!="):
                rule = FilterRule(
                    column=column, operator=operator,
                    value=str(self.input_num.value()),
                    is_numeric=True
                )
            else:
                text = self.input_text.text().strip()
                if not text:
                    QMessageBox.information(self, "Value Required", "Please enter a value for the filter.")
                    return
                rule = FilterRule(column=column, operator=operator, value=text, is_numeric=False)

            self.filter_manager.add_filter(rule)
            self._update_filters_display()
            self._update_values_list()

        except Exception as e:
            QMessageBox.warning(self, "Filter Error", f"Could not add filter: {e}")

    def _add_value_filter(self):
        """Add filter from selected value in values list."""
        item = self.values_list.currentItem()
        if not item:
            QMessageBox.information(self, "Select Value", "Please select a value from the list first.")
            return

        value = item.data(Qt.UserRole) or ""
        column = self.column_combo.currentText()
        is_numeric = column in NUMERIC_COLUMNS

        rule = FilterRule(column=column, operator="=", value=value, is_numeric=is_numeric)
        self.filter_manager.add_filter(rule)
        self._update_filters_display()
        self._update_values_list()

    def _reset_filters(self):
        """Clear all filters."""
        self.filter_manager.clear_filters()
        self._update_filters_display()
        self._update_values_list()

        if self.auto_filter_button.isChecked():
            self.auto_filter_button.setChecked(False)
            self.auto_filter_button.setText("Auto-Filter: OFF")

        # Also clear HPO filters
        if hasattr(self, 'hpo_filter_widget'):
            self.hpo_filter_widget._clear_hpo_filters()

    def _update_filters_display(self):
        """Update the display of active filters."""
        # Clear existing filter widgets
        while self.filters_vbox.count() > 1:
            item = self.filters_vbox.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

        # Add current filters
        for rule in self.filter_manager.filters:
            filter_widget = self._create_filter_widget(rule)
            self.filters_vbox.insertWidget(self.filters_vbox.count() - 1, filter_widget)

    def _create_filter_widget(self, rule: FilterRule) -> QWidget:
        """Create widget to display a single filter rule."""
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)

        # Rule description and remove button
        row_widget = QWidget()
        row_layout = QHBoxLayout(row_widget)
        row_layout.setContentsMargins(0, 0, 0, 0)

        description_label = QLabel(rule.describe())
        remove_button = QPushButton("Remove")
        remove_button.setFixedWidth(90)
        remove_button.clicked.connect(lambda: self._remove_filter(rule))

        row_layout.addWidget(description_label, 1)
        row_layout.addWidget(remove_button, 0)

        # Separator line
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)

        layout.addWidget(row_widget)
        layout.addWidget(separator)

        return container

    def _remove_filter(self, rule: FilterRule):
        """Remove a specific filter rule."""
        self.filter_manager.remove_filter(rule)
        self._update_filters_display()
        self._update_values_list()

    def _update_values_list(self):
        """Update the list of distinct values for the selected column."""
        self.values_list.clear()
        column = self.column_combo.currentText()

        if not column or column not in ALL_COLUMNS:
            return

        col_idx = ALL_COLUMNS.index(column)
        value_counts = {}

        # Count visible values
        for row in range(self.table.rowCount()):
            if self.table.isRowHidden(row):
                continue

            item = self.table.item(row, col_idx)
            value = item.text() if item else ""
            value_counts[value] = value_counts.get(value, 0) + 1

        # Sort values appropriately
        if column in NUMERIC_COLUMNS:
            values = sorted(value_counts.keys(), key=lambda x: self._safe_numeric_sort(x))
        else:
            values = sorted(value_counts.keys(), key=str.lower)

        # Add to list widget
        for value in values:
            count = value_counts[value]
            display_text = f"{value}  ({count})" if value else f"<empty>  ({count})"
            item = QListWidgetItem(display_text)
            item.setData(Qt.UserRole, value)
            self.values_list.addItem(item)

    def _safe_numeric_sort(self, value: str) -> float:
        """Safely convert value to number for sorting."""
        try:
            return float(value) if value.strip() else float('inf')
        except ValueError:
            return float('inf')

    def _filter_values_list(self):
        """Filter the values list based on search text."""
        search_text = self.values_search.text().lower()

        for i in range(self.values_list.count()):
            item = self.values_list.item(i)
            value = item.data(Qt.UserRole) or ""
            display_value = value if value else "<empty>"
            item.setHidden(search_text and search_text not in display_value.lower())

    # Column visibility management
    def _build_column_menu(self):
        """Build the column visibility menu."""
        self.columns_menu.clear()

        # Show all checkbox
        show_all_checkbox = QCheckBox("Show all", self)
        show_all_checkbox.setChecked(all(self.column_visibility.values()))
        show_all_checkbox.toggled.connect(self._toggle_all_columns)

        show_all_action = QWidgetAction(self.columns_menu)
        show_all_action.setDefaultWidget(show_all_checkbox)
        self.columns_menu.addAction(show_all_action)
        self.columns_menu.addSeparator()

        # Individual column checkboxes
        self.column_actions = {}
        for column in ALL_COLUMNS:
            checkbox = QCheckBox(column, self)
            checkbox.setChecked(self.column_visibility.get(column, True))
            checkbox.toggled.connect(lambda checked, col=column: self._toggle_column(col, checked))

            action = QWidgetAction(self.columns_menu)
            action.setDefaultWidget(checkbox)
            self.columns_menu.addAction(action)
            self.column_actions[column] = checkbox

    def _toggle_all_columns(self, checked: bool):
        """Toggle visibility of all columns."""
        for column in ALL_COLUMNS:
            self.column_visibility[column] = checked
            if column in self.column_actions:
                self.column_actions[column].blockSignals(True)
                self.column_actions[column].setChecked(checked)
                self.column_actions[column].blockSignals(False)
        self._apply_column_visibility()

    def _toggle_column(self, column: str, checked: bool):
        """Toggle visibility of a specific column."""
        self.column_visibility[column] = checked
        self._apply_column_visibility()

    def _apply_column_visibility(self):
        """Apply column visibility settings to the table."""
        for i, column in enumerate(ALL_COLUMNS):
            self.table.setColumnHidden(i, not self.column_visibility.get(column, True))

    # Table population and data management
    def _populate_table(self):
        """Populate table with VCF data."""
        self.table.setSortingEnabled(False)
        self.table.clearContents()
        self.table.setRowCount(0)

        if not self.variants:
            self.table.setSortingEnabled(True)
            return

        # Process dynamic FORMAT columns
        self._process_format_columns()

        # Update table structure
        self.table.setColumnCount(len(ALL_COLUMNS))
        self.table.setHorizontalHeaderLabels(ALL_COLUMNS)
        self.table.setRowCount(len(self.variants))

        # Populate data
        col_mapping = {name: idx for idx, name in enumerate(self.headers)}
        format_idx = col_mapping.get("FORMAT")
        sample_idx = format_idx + 1 if format_idx is not None and format_idx + 1 < len(self.headers) else None

        for row_idx, variant_fields in enumerate(self.variants):
            self._populate_table_row(row_idx, variant_fields, col_mapping, format_idx, sample_idx)

        # Update UI components
        self._refresh_ui_after_load()

    def _process_format_columns(self):
        """Process FORMAT fields to create dynamic columns."""
        if not self.variants:
            return

        col_mapping = {name: idx for idx, name in enumerate(self.headers)}
        format_idx = col_mapping.get("FORMAT")

        if format_idx is None:
            return

        # Collect all FORMAT keys
        format_keys = []
        seen_keys = set()

        for variant_fields in self.variants:
            if format_idx < len(variant_fields) and variant_fields[format_idx]:
                for key in variant_fields[format_idx].split(":"):
                    if key and key not in seen_keys:
                        seen_keys.add(key)
                        format_keys.append(key)

        # Add new columns to globals - INSERT before annotation columns instead of appending
        self.dynamic_format_columns = format_keys

        # Find the insertion point (before annotation columns)
        insertion_point = len(BASE_COLUMNS)

        for key in format_keys:
            if key not in ALL_COLUMNS:
                ALL_COLUMNS.insert(insertion_point, key)  # Insert instead of append
                insertion_point += 1  # Increment for next insertion
                self.column_visibility[key] = (key in ["GT", "AF"])
                # Add numeric columns heuristically
                if key in {"DP", "AO", "GQ", "AF"}:
                    NUMERIC_COLUMNS.add(key)

    def _populate_table_row(self, row_idx: int, variant_fields: List[str],
                            col_mapping: Dict[str, int], format_idx: Optional[int],
                            sample_idx: Optional[int]):
        """Populate a single table row."""
        # Base VCF columns
        for col_name in BASE_COLUMNS:
            value = self._get_column_value(col_name, variant_fields, col_mapping, format_idx, sample_idx)
            col_idx = ALL_COLUMNS.index(col_name)
            self.table.setItem(row_idx, col_idx, QTableWidgetItem(value))

        # Dynamic FORMAT columns
        if format_idx is not None and sample_idx is not None:
            format_mapping = self._parse_format_values(variant_fields, format_idx, sample_idx)
            for key in self.dynamic_format_columns:
                value = format_mapping.get(key, "")
                col_idx = ALL_COLUMNS.index(key)
                self.table.setItem(row_idx, col_idx, QTableWidgetItem(value))

        # Empty annotation columns
        for col_name in ANNOTATION_COLUMNS:
            col_idx = ALL_COLUMNS.index(col_name)
            self.table.setItem(row_idx, col_idx, QTableWidgetItem(""))

        # Empty HPO columns
        for col_name in HPO_COLUMNS:
            col_idx = ALL_COLUMNS.index(col_name)
            self.table.setItem(row_idx, col_idx, QTableWidgetItem(""))

    def _get_column_value(self, col_name: str, variant_fields: List[str],
                          col_mapping: Dict[str, int], format_idx: Optional[int],
                          sample_idx: Optional[int]) -> str:
        """Get value for a specific column."""
        if col_name == "Selector Values":
            return variant_fields[sample_idx] if sample_idx and sample_idx < len(variant_fields) else ""
        elif col_name == "FORMAT":
            return variant_fields[format_idx] if format_idx and format_idx < len(variant_fields) else ""
        elif col_name == "#CHROM":
            idx = col_mapping.get("CHROM", col_mapping.get("#CHROM"))
        else:
            idx = col_mapping.get(col_name)

        return variant_fields[idx] if idx is not None and idx < len(variant_fields) else ""

    def _parse_format_values(self, variant_fields: List[str], format_idx: int,
                             sample_idx: int) -> Dict[str, str]:
        """Parse FORMAT field values."""
        format_keys = variant_fields[format_idx].split(":") if format_idx < len(variant_fields) else []
        format_values = variant_fields[sample_idx].split(":") if sample_idx < len(variant_fields) else []

        # Pad values to match keys
        while len(format_values) < len(format_keys):
            format_values.append("")

        return dict(zip(format_keys, format_values))

    def _populate_hpo_data(self, row_idx: int, gene_symbol: str):
        """Populate HPO-related data for a gene."""
        try:
            # Get HPO terms for this gene
            hpo_terms = self.hpo_manager.get_hpo_terms_for_gene(gene_symbol)

            if hpo_terms:
                # Get details for HPO terms
                term_names = []
                term_descriptions = []

                for hpo_id in list(hpo_terms)[:5]:  # Limit to first 5 terms
                    term_data = self.hpo_manager.get_hpo_term_details(hpo_id)
                    if term_data:
                        term_names.append(f"{hpo_id}:{term_data.get('name', '')}")
                        term_descriptions.append(term_data.get('definition', '')[:50] + '...' if len(term_data.get('definition', '')) > 50 else term_data.get('definition', ''))

                # Calculate a simple score based on number of associated HPO terms
                hpo_score = min(len(hpo_terms) / 10.0, 1.0)  # Normalize to 0-1 range

                # Update table cells
                hpo_terms_col = ALL_COLUMNS.index("HPO Terms")
                hpo_desc_col = ALL_COLUMNS.index("HPO Descriptions")
                hpo_score_col = ALL_COLUMNS.index("Gene-HPO Score")

                self.table.setItem(row_idx, hpo_terms_col, QTableWidgetItem("; ".join(term_names)))
                self.table.setItem(row_idx, hpo_desc_col, QTableWidgetItem("; ".join(term_descriptions)))
                self.table.setItem(row_idx, hpo_score_col, QTableWidgetItem(f"{hpo_score:.3f}"))
            else:
                # Set empty values
                hpo_terms_col = ALL_COLUMNS.index("HPO Terms")
                hpo_desc_col = ALL_COLUMNS.index("HPO Descriptions")
                hpo_score_col = ALL_COLUMNS.index("Gene-HPO Score")

                self.table.setItem(row_idx, hpo_terms_col, QTableWidgetItem(""))
                self.table.setItem(row_idx, hpo_desc_col, QTableWidgetItem(""))
                self.table.setItem(row_idx, hpo_score_col, QTableWidgetItem(""))

        except Exception as e:
            print(f"Error populating HPO data for {gene_symbol}: {e}")

    def _refresh_ui_after_load(self):
        """Refresh UI components after loading data."""
        self.table.setSortingEnabled(True)
        self._build_column_menu()
        self.column_combo.clear()
        self.column_combo.addItems(ALL_COLUMNS)
        self._update_operator_options()
        self._update_values_list()
        self._apply_column_visibility()

    def _update_header_inspector(self, file_path: str):
        """Update header inspector with VCF header information."""
        try:
            header_lines = VCFReader.collect_header_lines(file_path)
            self.header_inspector.set_header_lines(header_lines)
        except Exception as e:
            print(f"Could not load header information: {e}")


    def _load_hpo_files(self):
        """Load HPO files for local filtering."""
        # Get hp.obo file
        hp_obo_path, _ = QFileDialog.getOpenFileName(
            self, "Select hp.obo file", "", "OBO files (*.obo);;All files (*)"
        )
        if not hp_obo_path:
            return

        # Get phenotype_to_genes.txt file
        genes_path, _ = QFileDialog.getOpenFileName(
            self, "Select phenotype_to_genes.txt file", "", "Text files (*.txt);;All files (*)"
        )
        if not genes_path:
            return

        # Show progress dialog
        progress_dialog = QMessageBox(self)
        progress_dialog.setWindowTitle("Loading HPO Files")
        progress_dialog.setText("Loading HPO files... This may take a moment.")
        progress_dialog.setStandardButtons(QMessageBox.NoButton)
        progress_dialog.setModal(True)  # Make it modal
        progress_dialog.show()
        QApplication.processEvents()

        try:
            # Load files
            self.hpo_manager.load_hpo_files(hp_obo_path, genes_path)

            progress_dialog.hide()  # Use hide() instead of close()
            progress_dialog.deleteLater()  # Ensure cleanup

            if self.hpo_manager.loaded:
                # Update HPO widget status
                if hasattr(self, 'hpo_filter_widget'):
                    self.hpo_filter_widget.set_hpo_loaded(True)

                QMessageBox.information(
                    self, "Success",
                    f"HPO files loaded successfully!\n\n"
                    f"Terms: {len(self.hpo_manager.hpo_terms)}\n"
                    f"Gene associations: {len(self.hpo_manager.gene_hpo_map)}"
                )
            else:
                QMessageBox.warning(self, "Error", "Failed to load HPO files.")

        except Exception as e:
            progress_dialog.hide()  # Use hide() instead of close()
            progress_dialog.deleteLater()  # Ensure cleanup
            QMessageBox.critical(self, "Error", f"Error loading HPO files:\n{str(e)}")

    def _debug_gene_column(self):
        """Debug method to see what's in the Gene column."""
        gene_col_idx = ALL_COLUMNS.index("Gene")
        genes_found = set()

        for row in range(min(20, self.table.rowCount())):  # Check first 20 rows
            gene_item = self.table.item(row, gene_col_idx)
            gene_value = gene_item.text().strip() if gene_item else ""
            if gene_value:
                genes_found.add(gene_value)

        print(f"Genes found in table: {list(genes_found)}")
        return genes_found

    def _populate_all_hpo_data_with_count(self):
        """Populate HPO data for all rows that have gene information and return count."""
        if not self.hpo_manager.loaded:
            return 0

        gene_col_idx = ALL_COLUMNS.index("Gene")
        populated_count = 0

        for row in range(self.table.rowCount()):
            gene_item = self.table.item(row, gene_col_idx)
            if gene_item and gene_item.text().strip():
                gene_symbol = gene_item.text().strip()
                self._populate_hpo_data(row, gene_symbol)
                populated_count += 1

        print(f"Populated HPO data for {populated_count} rows")
        return populated_count

    def _populate_all_hpo_data(self):
        """Populate HPO data for all rows that have gene information."""
        return self._populate_all_hpo_data_with_count()

    def _toggle_auto_filter(self):
        """Toggle auto-filter on/off."""
        if self.auto_filter_button.isChecked():
            self.auto_filter_button.setText("Auto-Filter: ON")
            self._apply_standard_filters()
        else:
            self.auto_filter_button.setText("Auto-Filter: OFF")
            self.filter_manager.clear_automatic_filters()
            self._update_filters_display()
            self._update_values_list()

    def _apply_standard_filters(self):
        """Apply standard quality filters automatically."""
        if not self.auto_filter_button.isChecked():
            return

        self.filter_manager.clear_automatic_filters()

        # Define standard filters
        standard_filters = []

        # 1. Allele Frequency >= 1% (0.01)
        # Check both possible frequency columns
        for freq_col in ["gnomAD AF", "Max Allele Freq"]:
            if freq_col in ALL_COLUMNS:
                freq_filter = FilterRule(
                    column=freq_col,
                    operator=">=",
                    value="0.01",
                    is_numeric=True,
                    is_automatic=True
                )
                standard_filters.append(freq_filter)
                break  # Only add one frequency filter

        # 2. Quality >= 20
        if "QUAL" in ALL_COLUMNS:
            qual_filter = FilterRule(
                column="QUAL",
                operator=">=",
                value="20",
                is_numeric=True,
                is_automatic=True
            )
            standard_filters.append(qual_filter)

        # 3. Impact at least "LOW" (exclude MODIFIER)
        if "Impact" in ALL_COLUMNS:
            impact_filter = FilterRule(
                column="Impact",
                operator="not contains",
                value="MODIFIER",
                is_numeric=False,
                is_automatic=True
            )
            standard_filters.append(impact_filter)

        # Apply all standard filters
        for rule in standard_filters:
            self.filter_manager.add_filter(rule)

        # Update UI
        self._update_filters_display()
        self._update_values_list()

        print(f"Applied {len(standard_filters)} standard filters")

#Reader of VCF logic
class VCFReader:
    """Utility class for reading VCF files."""

    @staticmethod
    def read_vcf(file_path: str) -> Tuple[List[str], List[List[str]]]:
        """Read VCF file and return headers and variant data."""
        import gzip

        opener = gzip.open if file_path.endswith(".gz") else open
        headers = []
        variants = []

        with opener(file_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()

                if line.startswith("##"):
                    continue
                elif line.startswith("#CHROM"):
                    headers = line.lstrip("#").split("\t")
                    continue

                fields = line.split("\t")
                if len(fields) >= 8:  # Minimum VCF fields
                    variants.append(fields)

        # Normalize chromosome column name
        if "CHROM" not in headers and "#CHROM" in headers:
            headers[headers.index("#CHROM")] = "CHROM"

        return headers, variants

    @staticmethod
    def collect_header_lines(file_path: str) -> List[str]:
        """Collect relevant header lines for inspector."""
        import gzip

        opener = gzip.open if file_path.endswith(".gz") else open
        header_lines = []

        with opener(file_path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()

                if line.startswith('##'):
                    if any(line.startswith(f'##{prefix}=') for prefix in ['INFO', 'FILTER', 'FORMAT']):
                        header_lines.append(line)
                elif line.startswith('#CHROM'):
                    break

        return header_lines


# Application entry point
def main():
    """Main application entry point."""
    app = QApplication(sys.argv)
    window = AppWindow()
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()

# Dependencies
# PyQt5==5.15.11
# PyQt5-sip==12.17.0
# PyQt5-Qt5==5.15.17
# requests>=2.25.0